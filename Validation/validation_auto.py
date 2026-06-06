import sys
import Sofa
import SofaRuntime
import meshio
import trimesh

import math
from scipy.spatial import ConvexHull
import numpy as np
import pandas as pd
import csv

import multiprocessing
from tqdm import tqdm 

import os

SofaRuntime.importPlugin("Sofa.Component")

mshFile = '../assets/gallbladder/lvl3/gb_3.msh'

# '../assets/gallbladder/lvl4/gb_4_v.msh' lvl 4

fixedPlane = [-1, 0, 0]

points = []
points_surface = []

vid_disp = [0.001775778,0.002124548,0.002047679,0.002211965,0.001409598,0.002723787,0.002192428,0.002657578,0.000621132,0.002104937,0.002806445,0.0017678,0.002998532,0.002593567,0.002911985,0.003785514,0.001904528,0.004912242,0.006542069,0.002205436,0.004879828,0.010454136,0.00939024,0.005144637,0.008896377,0.008690593,0.005363989,0.006787468,0.005223584,0.006024561,0.004115292,0.004706141,0.001899449,0.004629839,0.004569925,0.003842827,0.007358106,0.002337578,0.002390833,0.000717505,0.002357316
]


# [0.001185426, 0.001104278, 0.001190776, 0.001374173, 0.00181308, 0.001556994, 0.00127625, 0.001116709, 0.00163796, 0.004855458, 0.00579759, 0.005726778, 0.005546219, 0.003946628, 0.004952723, 0.005910661, 0.006671463, 0.007027979, 0.006852572, 0.007765475, 0.007926983, 0.007881653, 0.005870132, 0.00750954, 0.008525874, 0.008470684, 0.007926572, 0.003377693, 0.010325625, 0.009721619, 0.00795634, 0.007032534, 0.003285961, 0.006739244, 0.003080373]



mesh = meshio.read(mshFile)
points = mesh.points

def mesh_setup():
    #directly modifies global variables 
    global points_surface, cells_surface, mesh_center
    surface = set()
    for cell_block in mesh.cells:
        if cell_block.type in ['triangle','quad','line']:
            surface.update(cell_block.data.flatten()) 

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
        
    xBar = xMaxT - (xMaxT - xMinT)/2.5 #for back (liver-wall) springs 5
    xBarF = xMaxT - (xMaxT- xMinT)/7 #for fixed constraints 7
    zBarTop = zMaxT - (zMaxT - zMinT)/25 #for top springs
    zBarBot = zMinT + (zMaxT - zMinT)/10 #for bottom springs

    index = 0
    for point in points_surface: 
        if point[0] >= xBar or point[2] >= zBarTop or point[2] <= zBarBot:
            constrainedIndices.append(index)
        if point[0] >= xBarF:
            fixedIndices.append(index)
        index += 1

    return constrainedIndices, fixedIndices

def createScene(root, retract_ind, pull_vector, pull_force, points_surface, stiffness):
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
    ])

    root.addObject('DefaultVisualManagerLoop')

    root.addObject('CollisionPipeline',verbose='1')
    root.addObject('BruteForceBroadPhase')
    root.addObject('BVHNarrowPhase')
    root.addObject('CollisionResponse', name='collision response', response='PenalityContactForceField')
    root.addObject('DiscreteIntersection')

    root.addObject('DefaultAnimationLoop')

    root.addObject('EulerImplicitSolver', rayleighStiffness='0.1', rayleighMass='0.1')
    root.addObject('CGLinearSolver', iterations='30', tolerance="1e-8", threshold="1e-8")

    gallbladder = root.addChild('Gallbladder')
    gallbladder.addObject('MeshGmshLoader',name='volMesh',filename=mshFile,rotation='0 0 0')
    gallbladder.addObject('TetrahedronSetTopologyContainer',name='volTopo',src='@volMesh')
    gallbladder.addObject('MechanicalObject',name='dofs',src='@volMesh',template="Vec3d")
    gallbladder.addObject('UniformMass',totalMass='0.04')
    gallbladder.addObject('TetrahedronFEMForceField',name='FEM',youngModulus=stiffness,poissonRatio='0.4',method='large',computeVonMisesStress='1')
    gallbladder.addObject('TetrahedronSetGeometryAlgorithms')
    gallbladder.addObject('RestShapeSpringsForceField',name='conTissue',points=np.array(constrainedIndices),stiffness=12,drawSpring=1)
    # gallbladder.addObject('FixedProjectiveConstraint', indices=fixedIndices)
    
    gallbladder.addObject('SphereROI',name='forceROI',centers=points_surface[retract_ind], radii=0.0075,position='@dofs.position')
    indForceROI = gallbladder.forceROI.indices.value
    forceV = np.multiply(pull_force/len(indForceROI), pull_vector)
    forceVTot = [] 
    for i in indForceROI:
        forceVTot.extend(forceV)
    
    gallbladder.addObject('ConstantForceField',name='pull',indices=indForceROI,forces=[float(x) for x in forceVTot])


def run_simulation(steps, retract_ind, pull_vector, pull_force, points_surface, stiffness):
    root = Sofa.Core.Node('root')
    createScene(root, retract_ind, pull_vector, pull_force, points_surface, stiffness) #now includes generateConstraints
    Sofa.Simulation.init(root)
    Sofa.Simulation.animateNSteps(root, steps, 0.02)
    return root

def retractions_at_point(retract_ind, pull_force, pull_vector, points_surface, stiffness, proj_pts):
    root = run_simulation(250, retract_ind, pull_vector, pull_force, points_surface, stiffness)
    sofa_disp = []
    for i in range(len(proj_pts)):
        if proj_pts[i] == None:
            sofa_disp.append(vid_disp[i]) 
            continue
        else:
            o_pos = np.array(points[proj_pts[i]])
            n_pos = np.array(root.Gallbladder.dofs.position.value[proj_pts[i]])
        disp_m = np.linalg.norm(n_pos - o_pos)
        sofa_disp.append(disp_m)
    raw_error = []
    abs_error = []
    abs_error_pcnt = []
    for ind in range(len(sofa_disp)):
        raw_error.append(sofa_disp[ind]-vid_disp[ind])
        abs_error.append(abs(sofa_disp[ind]-vid_disp[ind]))
        abs_error_pcnt.append(abs_error[ind]/vid_disp[ind])
    raw_error_avg = sum(raw_error)/len(raw_error)
    abs_error_avg = sum(abs_error)/len(abs_error)
    abs_error_pcnt_avg = sum(abs_error_pcnt)/len(abs_error_pcnt)
    # relative_error_pct =abs_error/sofa_disp


    results_df = pd.DataFrame({
        'vid_disp': vid_disp,
        'sofa_disp': sofa_disp,
        'raw_error': raw_error,
        'abs_error': abs_error,
        'abs_error_pcnt': abs_error_pcnt
    })

    results_df.to_csv(f'{stiffness}Pa Retraction Results.csv')

    Sofa.Simulation.unload(root)
    del root
    return raw_error_avg, abs_error_avg, abs_error_pcnt_avg


def main():
    mesh_setup()
    
    
    retract_ind = 515
    # lvl 4 retractionindex: 952
    #lvl 1 retraction index: 559
    force = 7
    dir_v = [0, 0.9806, -0.1961]
    proj_pts = [1302, 855, 694, 1341, 876, 1312, 700, 1122, 1365, 623, 1400, 1256, 947, 514, 1053, 803, 66, 1086, 194, 73, 633, 938, 201, 272, 1178, 509, 1301, 1307, 408, 859, 310, 3, 198, 430, 526, 569, 328, 20, 663, 682, 775]
    last_abs_error_pcnt_avg = 0

    low = 10000
    high = 100000
    mid = 0
    prev_abs_error = 0
    # for iter in range(8):
    #     mid = (high+low)//2
    #     if iter == 0:
    #         raw_error_avg, abs_error_avg, abs_error_pcnt_avg = retractions_at_point(retract_ind, force, dir_v, points_surface, mid,proj_pts)
    #         stiff = mid
    #     else:
    #         stiff_L = (low+mid)//2
    #         reaL, aeaL, aepaL = retractions_at_point(retract_ind, force, dir_v, points_surface, stiff_L,proj_pts)

    #         stiff_H = (mid+high)//2
    #         reaH, aeaH, aepaH = retractions_at_point(retract_ind, force, dir_v, points_surface, stiff_H,proj_pts)

    #         if aepaL < aepaH: 
    #             high = mid
    #             raw_error_avg, abs_error_avg, abs_error_pcnt_avg = reaL, aeaL, aepaL
    #             stiff = stiff_L
    #         else: 
    #             low = mid
    #             raw_error_avg, abs_error_avg, abs_error_pcnt_avg = reaH, aeaH, aepaH
    #             stiff = stiff_H
        
    #     print(stiff,raw_error_avg, abs_error_avg, abs_error_pcnt_avg)
            

    for stiffness in np.linspace(115000, 125000, 12):
        raw_error_avg, abs_error_avg, abs_error_pcnt_avg, = retractions_at_point(retract_ind, force, dir_v, points_surface, int(stiffness), proj_pts)
        print(stiffness,raw_error_avg, abs_error_avg, abs_error_pcnt_avg)
    #     if iter == 0:
    #         if raw_error_avg >= 0:
    #             high = mid
    #             stiff_incr = False
    #         else:
    #             low = mid
    #     else:
    #         if abs_error_pcnt_avg > last_abs_error_pcnt_avg:
    #             stiff_incr = not(stiff_incr)
    #         if stiff_incr:
    #             low = mid
    #         else:
    #             high = mid
    #     last_abs_error_pcnt_avg = abs_error_pcnt_avg
    
    # print('failsafe')
    # low = 10000
    # high = 100000
    # for iter in range(8):
    #     mid = (low + high)//2
    #     print(mid)
    #     raw_error_avg, abs_error_avg, abs_error_pcnt_avg = retractions_at_point(retract_ind, force, dir_v, points_surface, mid,proj_pts)
    #     print(raw_error_avg, abs_error_avg, abs_error_pcnt_avg)
        
    #     if iter == 0:
    #         if raw_error_avg <= 0:
    #             high = mid
    #             stiff_incr = False
    #         else:
    #             low = mid
    #     else:
    #         if abs_error_pcnt_avg > last_abs_error_pcnt_avg:
    #             stiff_incr = not(stiff_incr)
    #         if stiff_incr:
    #             low = mid
    #         else:
    #             high = mid
    #     last_abs_error_pcnt_avg = abs_error_pcnt_avg

    # print("raw")
    # low = 10000
    # high = 100000
    # for iter in range(8):
    #     mid = (low + high)//2
    #     print(mid)
    #     raw_error_avg, abs_error_avg, abs_error_pcnt_avg = retractions_at_point(retract_ind, force, dir_v, points_surface, mid,proj_pts)
    #     print(raw_error_avg, abs_error_avg, abs_error_pcnt_avg)
        
    #     if raw_error_avg <= 0:
    #         high = mid
    #     else:
    #         low = mid

if __name__ == '__main__':
    main()