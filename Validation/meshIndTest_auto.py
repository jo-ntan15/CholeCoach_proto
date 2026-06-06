import sys
import Sofa
import SofaRuntime
import meshio
import trimesh

import math
import numpy as np
import pandas as pd

import multiprocessing

SofaRuntime.importPlugin("Sofa.Component")

#visual check block
import Sofa.Gui
SofaRuntime.importPlugin("Sofa.GL.Component")

mesh_files = [
    # '../assets/gallbladder/lvl4-2/gb_4-2_s2.msh',
    'meshIndTest_meshes/lvl4-2/gb_4-2_s1.msh',
    'meshIndTest_meshes/lvl4-2/gb_4-2_s.msh',
    'meshIndTest_meshes/lvl4-2/gb_level4_2.msh',
    'meshIndTest_meshes/lvl4-2/gb_4-2_v.msh',
    'meshIndTest_meshes/lvl4-2/gb_4-2_v1.msh',
    'meshIndTest_meshes/lvl4-2/gb_4-2_v2.msh',
    'meshIndTest_meshes/lvl4-2/gb_4-2_v3.msh'
]

fixedPlane = [-1, 0, 0]

stiffness = 50000

points = []
points_surface = []
points_ref = meshio.read(mesh_files[2]).points
cells_surface = []


def mesh_setup(file):
    #directly modifies global variables 
    global points, points_surface, cells_surface
    mesh = meshio.read(file)
    points = mesh.points
    surface = set()
    for cell_block in mesh.cells:
        if cell_block.type in ['triangle','quad','line']:
            surface.update(cell_block.data.flatten()) 
            #cell_block.data contains all connections, surface (set) will include all points that are a part of a surface cell
            cell_type = cell_block.type
            connections = cell_block.data
            for cell in connections:
                cells_surface.append([cell_type, list(cell)])

    points_surface = mesh.points[list(surface)]

def createPlane(vLst):
    vNp = []
    for v in vLst:
        vNp.append(np.array(v))
    v1 = vNp[2] - vNp[0]
    v2 = vNp[1] - vNp[0]
    normal = np.cross(v1, v2)
    D = -np.dot(normal, vNp[0])
    return (*normal, D) # * unpacks normal into three independent floats

def generateConstraints(points_surface):
    #points_surface (array) -> constraintIndices, fixedIndices (arrays)
    constrainedIndices = []
    fixedIndices = []

    xMinT = np.min(points_surface[:,0])
    xMaxT = np.max(points_surface[:,0])
    zMinT = np.min(points_surface[:,2])
    zMaxT = np.max(points_surface[:,2])
        
    xBar = xMaxT - (xMaxT - xMinT)/3 #for back (liver-wall) springs
    xBarF = xMaxT - (xMaxT- xMinT)/5 #for fixed constraints
    zBarTop = zMaxT - (zMaxT - zMinT)/25 #for top springs
    zBarBot = zMinT + (zMaxT - zMinT)/25 #for bottom springs

    index = 0
    for point in points_surface: 
        if point[0] >= xBar or point[2] >= zBarTop or point[2] <= zBarBot:
            constrainedIndices.append(index)
        if point[0] >= xBarF:
            fixedIndices.append(index)
        index += 1

    return constrainedIndices, fixedIndices

def createScene(root, mesh_file, points_surface):
    print(mesh_file)
    constrainedIndices, fixedIndices = generateConstraints(points_surface)
    # print('creating scene')
    root.dt = 0.02
    root.gravity = [0, -0.00001, 0]
    root.bbox = [[-0.04,-0.04,-0.04],[0.04,0.04,0.04]]
    
    root.addObject("RequiredPlugin", pluginName=[
        'Sofa.Component.Collision.Detection.Algorithm',
        'Sofa.Component.Collision.Detection.Intersection',
        'Sofa.Component.Collision.Geometry',
        'Sofa.Component.Collision.Response.Contact',
        'Sofa.Component.Constraint.Projective',
        'Sofa.Component.IO.Mesh',
        'Sofa.Component.LinearSolver.Iterative',
        'Sofa.Component.Mapping.Linear',
        'Sofa.Component.Mass',
        'Sofa.Component.ODESolver.Backward',
        'Sofa.Component.SolidMechanics.FEM.Elastic',
        'Sofa.Component.SolidMechanics.Spring',
        'Sofa.Component.StateContainer',
        'Sofa.Component.Topology.Container.Dynamic',
        'Sofa.Component.Visual',
        'Sofa.GL.Component.Rendering3D',
        'Sofa.Component.SceneUtility',
        'Sofa.GUI.Component',
        'Sofa.Component.Setting'
    ])

    root.addObject('DefaultVisualManagerLoop')

    root.addObject('CollisionPipeline',verbose='1')
    root.addObject('BruteForceBroadPhase')
    root.addObject('BVHNarrowPhase')
    root.addObject('CollisionResponse', name='collision response', response='PenalityContactForceField')
    root.addObject('DiscreteIntersection')

    root.addObject('DefaultAnimationLoop')

    root.addObject('EulerImplicitSolver', rayleighStiffness='0.1', rayleighMass='0.1')
    root.addObject('CGLinearSolver', iterations='25', tolerance="1e-7", threshold="1e-7")

    gallbladder = root.addChild('Gallbladder')
    gallbladder.addObject('MeshGmshLoader',name='volMesh',filename=mesh_file,rotation='0 0 0')
    gallbladder.addObject('TetrahedronSetTopologyContainer',name='volTopo',src='@volMesh')
    gallbladder.addObject('MechanicalObject',name='dofs',src='@volMesh',template="Vec3d")
    gallbladder.addObject('UniformMass',totalMass='0.04')
    gallbladder.addObject('TetrahedronFEMForceField',name='FEM',youngModulus=stiffness,poissonRatio='0.4',method='large',computeVonMisesStress='1')
    gallbladder.addObject('TetrahedronSetGeometryAlgorithms')
    gallbladder.addObject('RestShapeSpringsForceField',name='conTissue',points=np.array(constrainedIndices),stiffness=12,drawSpring=1)
    gallbladder.addObject('FixedProjectiveConstraint', indices=fixedIndices)
    
    gallbladder.addObject('SphereROI',name='stressROI1', centers=points_ref[337],radii=0.005,position='@dofs.position',drawROI = 1)
    gallbladder.addObject('SphereROI',name='stressROI2', centers=points_ref[1133],radii=0.005,position='@dofs.position', drawROI =1 )
    gallbladder.addObject('SphereROI',name='stressROI3', centers=points_ref[503],radii=0.005,position='@dofs.position', drawROI=1)
    indStressROI1 = gallbladder.stressROI1.indices.value
    indStressROI2 = gallbladder.stressROI2.indices.value
    indStressROI3 = gallbladder.stressROI3.indices.value

    gallbladder.addObject('SphereROI',name='forceROI',centers=points_ref[503], radii=0.005,position='@dofs.position')
    indForceROI = gallbladder.forceROI.indices.value
    
    forceV = np.multiply(5/len(indForceROI), [0,1,0])
    forceVTot = [] 
    for i in indForceROI:
        forceVTot.extend(forceV)
    
    gallbladder.addObject('ConstantForceField',name='pull',indices=indForceROI,forces=[float(x) for x in forceVTot])

    return indStressROI1, indStressROI2, indStressROI3

def run_simulation(steps, mesh_file, points_surface):
    root = Sofa.Core.Node('root')
    stressIndList1, stressIndList2, stressIndList3 = createScene(root, mesh_file, points_surface) #now includes generateConstraints
    Sofa.Simulation.init(root)
    Sofa.Simulation.animateNSteps(root, steps, 0.02)
    return root, stressIndList1, stressIndList2, stressIndList3

def chunkify_list(arr, n):
    return [arr[i::n] for i in range(n)]

def main():
    
    stress_wrt_mesh = {
        'mesh_size': [],
        'stress_avg': []
    }
    for gb in mesh_files:
        mesh_setup(gb)
        root, stressIndList1, stressIndList2, stressIndList3 = run_simulation(250, gb, points_surface)
        stress_avg1 = np.mean([root.Gallbladder.FEM.vonMisesPerNode[i] for i in stressIndList1])
        stress_avg2 = np.mean([root.Gallbladder.FEM.vonMisesPerNode[i] for i in stressIndList2])
        stress_avg3 = np.mean([root.Gallbladder.FEM.vonMisesPerNode[i] for i in stressIndList3])
        stress_avg = round((stress_avg1+stress_avg2+stress_avg3)/3,2)
        stress_wrt_mesh['mesh_size'].append(round(len(points)/54.6,2))
        stress_wrt_mesh['stress_avg'].append(stress_avg)

    results_df = pd.DataFrame(stress_wrt_mesh)
    results_df.to_csv('mesh_ind_test_4-2_auto.csv')

    # Sofa.Gui.GUIManager.Init("myscene", "qt")
    # Sofa.Gui.GUIManager.createGUI(root, __file__)

    # Sofa.Gui.GUIManager.SetDimension(2000, 1500)
    
    # Sofa.Gui.GUIManager.MainLoop(root)

if __name__ == '__main__':
    main()



