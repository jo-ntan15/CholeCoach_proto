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

mshFile = '../assets/gallbladder/lvl1/gb_1m.msh'

youngMod = 50000

points = []
points_surface = []
cells = []
cells_surface = []

taskPosLock = np.array([
    [-0.118058,-0.0467191,0.0532513],
    [-0.118058,-0.0467191,0.0532513],
    [-0.118058,-0.0467191,0.0532513],
])

taskLookLock = np.array([
    [2.52167,1.33669,-0.0701377],
    [2.52167,1.33669,-0.0701377],
    [2.52167,1.33669,-0.0701377],
])

ref_v = np.array([-0.0481, 0.9408, -0.3354])


mesh = meshio.read(mshFile)
points = mesh.points

def mesh_setup():
    #directly modifies global variables 
    global points_surface, cells_surface, mesh_center
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

def find_cells(index_interest, cells_surface):
    #find all cells containing a certain index of interest
    # -> list of indices of said cells within cells_surface

    cells_interest = []
    for cell_ind in range(len(cells_surface)):
        if index_interest in cells_surface[cell_ind][1]:
            cells_interest.append(cell_ind)
    
    return cells_interest

def normal_per_cell(cell, points_surface):
    #returns normal vector of cell
    #list of points making up said cell -> v_normal
    v1 = np.array(points_surface[cell[0]])
    v2 = np.array(points_surface[cell[1]])
    v3 = np.array(points_surface[cell[2]])
    v_normal = np.linalg.cross(v2-v1, v3-v2)

    return v_normal/(np.linalg.norm(v_normal))

def normal_at_point(index_interest, points_surface, cells_surface):
    # goes through every cell around point, calculates cell normals, takes average 
    # index of point of interest (retraction site) -> normal at point 
    cells_ind_list = find_cells(index_interest, cells_surface)
    v_normal_tot = np.array([0.0,0.0,0.0])
    for ind in cells_ind_list:
        v_normal_tot += normal_per_cell(cells_surface[ind][1], points_surface)
        
    return (v_normal_tot/np.linalg.norm(v_normal_tot))

def generate_unit_v(spacing):
    #given "spacing", generate list of unit vectors 
    #azimuthal and polar angles
    vectors = []
    for theta_deg in range(0, 181, spacing):
        theta = np.radians(theta_deg)
        if theta_deg in {0, 180}:
            # Add only one vector for the poles
            z = 1.0 if theta_deg == 0 else -1.0
            vectors.append((0.0, 0.0, z))
        else:
            for phi_deg in range(0, 360, spacing):
                phi = np.radians(phi_deg)
                x = np.sin(theta) * np.cos(phi)
                y = np.sin(theta) * np.sin(phi)
                z = np.cos(theta)
                vectors.append((x, y, z))
    return np.array(vectors)

def find_angle(v, ref):
    #takes np unit vector "v" and reference "ref"
    # Compute dot product and magnitudes
    dot_product = np.dot(v, ref)
    norm_u = np.linalg.norm(ref)
    norm_v = np.linalg.norm(v)
    
    # Avoid division by zero (if either vector is zero)
    if norm_u == 0 or norm_v == 0:
        raise ValueError("One of the vectors has zero magnitude.")
    
    # Clip to handle floating-point errors (ensure -1 ≤ cosθ ≤ 1)
    cos_theta = np.clip(dot_product / (norm_u * norm_v), -1.0, 1.0)
    theta_rad = np.arccos(cos_theta)
    theta_deg = np.degrees(theta_rad)
    
    return theta_deg

def filter_valid_angles(curr_retract_ind, dissec_target, normal, unit_v, bound, points_surface):
    # filter angles give point of retraction, dissection, normal vector of point of retraction

    curr_retract_pos = points_surface[curr_retract_ind]
    dissec_target_pos = points_surface[dissec_target]
    filter_plane_normal = dissec_target_pos - curr_retract_pos
    filter_plane_normal /= np.linalg.norm(filter_plane_normal)

    valid_v = []
    for v in unit_v:
        if find_angle(v, normal) < bound and find_angle(v, filter_plane_normal) > 90 and v[0] <= 0:
            valid_v.append(list(v))

    return valid_v
    
def get_valid_points(dissec_target, empty_radius, height_bound):
    # returns all eligible points wrt. the dissection target index
    # constrained by empty radius around target and height away from target (absolute z)
    valid_retract = []
    xMax = np.max(points_surface[:,0])
    xMin = np.min(points_surface[:,0])
    xBar = (xMax + xMin)/2
    index = 0
    for point in points_surface:
        if np.linalg.norm(point - points_surface[dissec_target]) > empty_radius and point[0] < xBar and abs(point[2]-points_surface[dissec_target][2]) <= height_bound:
            valid_retract.append((index, point))
        index += 1
    return valid_retract

def get_point_proj_area(root, camPos, camLookAt, points_surface, stressIndList):
    #get the projected area of ROI onto plane representing the viewport 
    #
    plane_norm = camLookAt - camPos
    plane_norm = plane_norm/np.linalg.norm(plane_norm)
    
    if not find_angle(np.array([1,0,0]),plane_norm):
        any_v = np.array([1,0,0])
    else:
        any_v = np.array([0,1,0])
    plane_x = np.cross(any_v, plane_norm)
    plane_x = plane_x/np.linalg.norm(plane_x)
    plane_y = np.cross(plane_x, plane_norm)
    plane_y = plane_y/np.linalg.norm(plane_y)
    points_pos_curr = root.Gallbladder.dofs.position.value

    proj_point_ls = []
    for point_ind in stressIndList:
        if point_ind >= len(points_surface): 
            continue
        point = points_pos_curr[point_ind]
        point_v = point - camPos

        new_point_x = np.dot(point_v,plane_x)
        new_point_y = np.dot(point_v,plane_y)
        proj_point_ls.append([new_point_x, new_point_y])

    proj_point_ls = np.array(proj_point_ls)
    hull = ConvexHull(proj_point_ls)
    area = hull.volume

    return area

def get_displacement(root, points_surface, indStressROI, ref_v):
    perp_disp_mag = []
    for index in indStressROI:
        if index >= len(points_surface): continue
        points_curr_pos = root.Gallbladder.dofs.position.value
        raw_disp = points_curr_pos[index] - points_surface[index]
        perp_disp_mag.append(np.linalg.norm(np.cross(raw_disp, ref_v)/np.linalg.norm(ref_v)))

    return sum(perp_disp_mag)/len(perp_disp_mag)

def createScene(root, dissec, valid_retract_pos, pull_vector, pull_force, points_surface):
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
    root.addObject('CGLinearSolver', iterations='30', tolerance="1e-7", threshold="1e-7")

    gallbladder = root.addChild('Gallbladder')
    gallbladder.addObject('MeshGmshLoader',name='volMesh',filename=mshFile,rotation='0 0 0')
    gallbladder.addObject('TetrahedronSetTopologyContainer',name='volTopo',src='@volMesh')
    gallbladder.addObject('MechanicalObject',name='dofs',src='@volMesh',template="Vec3d")
    gallbladder.addObject('UniformMass',totalMass='0.04')
    gallbladder.addObject('TetrahedronFEMForceField',name='FEM',youngModulus=youngMod,poissonRatio='0.4',method='large',computeVonMisesStress='1')
    gallbladder.addObject('TetrahedronSetGeometryAlgorithms')
    gallbladder.addObject('RestShapeSpringsForceField',name='conTissue',points=np.array(constrainedIndices),stiffness=12,drawSpring=1)
    gallbladder.addObject('FixedProjectiveConstraint', indices=fixedIndices)
    
    gallbladder.addObject('SphereROI',name='stressROI', centers=points_surface[dissec],radii=0.005,position='@dofs.position')
    indStressROI = gallbladder.stressROI.indices.value

    gallbladder.addObject('SphereROI',name='forceROI',centers=valid_retract_pos, radii=0.005,position='@dofs.position')
    indForceROI = gallbladder.forceROI.indices.value
    forceV = np.multiply(pull_force/len(indForceROI), pull_vector)
    forceVTot = [] 
    for i in indForceROI:
        forceVTot.extend(forceV)
    
    gallbladder.addObject('ConstantForceField',name='pull',indices=indForceROI,forces=[float(x) for x in forceVTot])

    return indStressROI

def run_simulation(steps, dissec, valid_retract_pos, pull_vector, pull_force, points_surface, cells_surface):
    root = Sofa.Core.Node('root')
    stressIndList = createScene(root, dissec, valid_retract_pos, pull_vector, pull_force, points_surface) #now includes generateConstraints
    Sofa.Simulation.init(root)
    Sofa.Simulation.animateNSteps(root, steps, 0.02)
    return root, stressIndList

def retractions_at_point(retract_ind, retract_pos, dissec_target, unit_v, pull_force, allowed_angle, points_surface, cells_surface, camPos, camLookAt, ref_v):
    results = []
    normal = normal_at_point(retract_ind, points_surface, cells_surface)
    valid_dir = filter_valid_angles(retract_ind, dissec_target, normal, unit_v, allowed_angle/2, points_surface)
    for pull_vector in valid_dir:
        root, stressIndList = run_simulation(200, dissec_target, retract_pos, pull_vector, pull_force, points_surface, cells_surface)
        stress_avg = np.mean([root.Gallbladder.FEM.vonMisesPerNode[i] for i in stressIndList])
        exposed_area = get_point_proj_area(root, camPos,camLookAt,points_surface,stressIndList)
        disp_avg = get_displacement(root, points_surface, stressIndList, ref_v)
        results.append({
            'retraction_ind': retract_ind,
            'pos_x': retract_pos[0],
            'pos_y': retract_pos[1],
            'pos_z': retract_pos[2],
            'dir_x': pull_vector[0],
            'dir_y': pull_vector[1],
            'dir_z': pull_vector[2],
            'area_exposed': exposed_area,
            'stress_induced': stress_avg,
            'displacement': disp_avg
        })
        Sofa.Simulation.unload(root)
        del root
    return results

def point_worker(args):
    return retractions_at_point(*args)

def main():
    mesh_setup()
    
    dissec_target = 596 #index of dissection goes here (640)
    activity_n = 0
    empty_radius = 0.01 #untouchable circular region around dissection site
    height_bound = 0.03 #height restriction away from dissection site    

    unit_vector_spacing = 30
    unit_vector_spacing_fine = 20
    allowed_angle = 150
    pull_force = 5

    print(points_surface[dissec_target])
    valid_retract = get_valid_points(dissec_target, empty_radius, height_bound)
    unit_v = generate_unit_v(unit_vector_spacing)
    fine_unit_v = generate_unit_v(unit_vector_spacing_fine)

    args_list = [
        (ind, pos, dissec_target, unit_v, pull_force, allowed_angle, points_surface, cells_surface, taskPosLock[activity_n], taskLookLock[activity_n], ref_v) 
        for ind, pos in valid_retract
    ]
    
    print(normal_at_point(64, points_surface, cells_surface))
    print((points_surface[64] - points_surface[596])/np.linalg.norm(points_surface[64] - points_surface[596]))
    print(len(valid_retract))

    with open('optimization potato.csv',mode='w',newline='') as csvfile:
        fieldnames = ['retraction_ind','pos_x','pos_y','pos_z','dir_x','dir_y','dir_z','area_exposed','stress_induced','displacement']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        with multiprocessing.get_context('spawn').Pool(4) as pool:
            for results in tqdm(pool.imap_unordered(point_worker, args_list, chunksize=1),total=len(args_list),desc='Simulating Retractions'):
                for res in results:
                    writer.writerow(res)
                csvfile.flush()

    results_df = pd.read_csv('optimization potato.csv')
    
    grouped = results_df.groupby('retraction_ind')
    max_indices = grouped['area_exposed'].idxmax()
    max_exposed_df = results_df.loc[max_indices]

    max_exposed_df.to_csv("optimization_results_coarse.csv")

    print('Coarse Ready')

    exposed_threshold = max_exposed_df['area_exposed'].quantile(0.80)
    fine_exposed_df = max_exposed_df[max_exposed_df['area_exposed'] >= exposed_threshold]
    
    fine_valid_retract = [
        (row.retraction_ind, np.array([row.pos_x, row.pos_y, row.pos_z]))
        for row in fine_exposed_df.itertuples(index=False)
    ]
    fine_args_list = [
        (ind, pos, dissec_target, fine_unit_v, pull_force, allowed_angle, points_surface, cells_surface, taskPosLock[activity_n], taskLookLock[activity_n], ref_v) 
        for ind, pos in fine_valid_retract
    ]

    
    with open('optimization potato.csv',mode='w',newline='') as csvfile:
        fieldnames = ['retraction_ind','pos_x','pos_y','pos_z','dir_x','dir_y','dir_z','area_exposed','stress_induced','displacement']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        with multiprocessing.get_context('spawn').Pool(4) as pool:
            for results in tqdm(pool.imap_unordered(point_worker, fine_args_list, chunksize=1),total=len(fine_args_list),desc='Simulating Retractions (fine)'):
                for res in results:
                    writer.writerow(res)
                csvfile.flush()
    
    fine_results_df = pd.read_csv('optimization potato.csv')
    
    fine_grouped = fine_results_df.groupby('retraction_ind')
    fine_max_indices = fine_grouped['area_exposed'].idxmax()
    fine_max_exposed_df = fine_results_df.loc[fine_max_indices]

    fine_max_exposed_df.to_csv("optimization_results_fine.csv")
        
    stlmesh = trimesh.load('../assets/gallbladder/lvl1/gb_1m.stl')

    stl_vert = stlmesh.vertices
    stl_faces = stlmesh.faces

    pd.DataFrame(stl_vert, columns=['x','y','z']).to_csv('stl_vertices.csv',index=False)
    pd.DataFrame(stl_faces,columns=['v1','v2','v3']).to_csv('stl_faces.csv',index=False)

if __name__ == '__main__':
    main()



