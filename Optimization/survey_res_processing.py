import numpy as np
import pandas as pd
import csv
import math
from scipy.spatial import KDTree
import meshio

mesh1 = meshio.read("../assets/gallbladder/lvl1/gb_1m.msh")
mesh2 = meshio.read("../assets/gallbladder/lvl2-2/gb_2-2m.msh")
mesh3 = meshio.read("../assets/gallbladder/lvl3/gb_3.msh")
mesh4 = meshio.read("../assets/gallbladder/lvl4/gb_4m.msh")
mesh5 = meshio.read("../assets/gallbladder/lvl4-2/gb_level4_2.msh")

def mesh_setup(mesh):
    #directly modifies global variables 
    global points_surface, cells_surface, mesh_center
    surface = set()
    for cell_block in mesh.cells:
        if cell_block.type in ['triangle','quad','line']:
            surface.update(cell_block.data.flatten()) 
    points_surface = mesh.points[list(surface)]
    return points_surface

    
vertices_ls = [mesh_setup(mesh1),mesh_setup(mesh2),mesh_setup(mesh3),mesh_setup(mesh4),mesh_setup(mesh5)]
results_df = pd.read_csv('survey_results_preprocess.csv')

def find_polar(dx, dy, dz):
    pol_rad = np.arcsin(dz)
    return np.degrees(pol_rad)

def find_azi(dx, dy, dz):
    pol_rad = np.arcsin(dz)
    azi_rad = np.arcsin(dx/np.cos(pol_rad)) 
    return np.degrees(azi_rad)

df_levels = {level: sub_df for level, sub_df in results_df.groupby("level")}

df_activities = []
for i in range(5):
    df_activities.append({activity: sub_df for activity, sub_df in df_levels[i+1].groupby("activity")})

for i in range(5):
    level = df_activities[i]
    for j in range(3):
        activity = level[j+1]
        print(activity)
        for retraction in activity.itertuples():
            pos = np.array([retraction.x, retraction.y, retraction.z])
            tree = KDTree(vertices_ls[i])
            dist, idx = tree.query(pos)

            new_pos = vertices_ls[i][idx]
            print(pos)
            print(new_pos)

            activity.at[retraction.Index, 'x'] = new_pos[0]
            activity.at[retraction.Index, 'y'] = new_pos[1]
            activity.at[retraction.Index, 'z'] = new_pos[2]

        activity.to_csv(f'retraction_results/L{i+1}A{j+1}.csv')

# df_interest = pd.read_csv('retraction_results/L2A2.csv')

# angles = []
# for row in df_interest.itertuples():
#     dir_arr = [row.dx, row.dy, row.dz]
#     azi = find_azi(*dir_arr)
#     polar = find_polar(*dir_arr)
#     angles.append([row.id, row.expert, azi, polar])

# angles_df = pd.DataFrame(
#     angles,
#     columns = ['id', 'expert', 'azi', 'polar']
# )

# print(angles_df)