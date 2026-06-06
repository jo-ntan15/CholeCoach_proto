import math
import numpy as np
import meshio

mshFile = '../assets/gallbladder/lvl4/gb_4_v.msh'
#change mshFile to the msh being validated (match it to the gb_survey_validation file)

mesh = meshio.read(mshFile)
points_msh = mesh.points

def angle(v1,v2):
    dot = np.dot(v1,v2)
    norm = math.sqrt(v1[0]**2 + v1[1]**2) * math.sqrt(v2[0]**2 + v2[1]**2)
    return math.acos(dot/norm)

def rotate_point(x, y, cx, cy, angle_rad):
    x -= cx
    y -= cy
    
    x_new = x * math.cos(angle_rad) - y * math.sin(angle_rad)
    y_new = x * math.sin(angle_rad) + y * math.cos(angle_rad)
    
    return (x_new + cx, y_new + cy)

points_given = [
    [-517,-19],   #-533 40 #508 17
    [-585,-473]    #-533 -436 #540 433
]
#replace points_given with the higherst and lowest points on the gallbladder
#this may need to be manually calibrated 

with open('validation_init.txt','r') as file:
    for line in file:
        points_given.append([-float(num) for num in line.strip().split()])
#paste the initial points here (do NOT include highest and lowest)

points_given = np.array(points_given)

zMin = 1
zMax = -1
index = zMaxI = zMinI = 0
for point in points_msh:
    zMin = min(zMin,point[2])
    zMax = max(zMax,point[2])
#identifying max min x y z indexes
for point in points_msh: 
    if point[2] == zMax:
        zMaxI = index
        print(index, 'yMinI update')
    if point[2] == zMin: 
        zMinI = index
        print(index, 'zMinI update')
    index += 1

points_template = np.array([points_msh[zMaxI][1:],points_msh[zMinI][1:]])

print(f'points given {points_given}')
print(f'template {points_template}')

A_A_diff = points_template[0] - points_given[0]
print(A_A_diff)
#translation
for point in points_given:
    point += A_A_diff

print(f'translated {points_given}')

#scale
given_diff = points_given[0]-points_given[1]
given_dist = math.sqrt(np.dot(given_diff,given_diff))
template_diff = points_template[0]-points_template[1]
template_dist = math.sqrt(np.dot(template_diff,template_diff))

scale_factor = template_dist/given_dist
print(scale_factor)
for ind in range(len(points_given)):
    if ind == 0: pass
    points_given[ind] = points_given[0] + (points_given[ind] - points_given[0])*scale_factor
print(f"scaled {points_given}")

rotation_angle = angle(points_given[0]-points_given[1],points_template[0]-points_template[1])

for ind in range(len(points_given)):
    if ind == 0: pass
    points_given[ind] = rotate_point(*points_given[ind],*points_given[0],rotation_angle) 
    # put negative sign in front of rotation_angle

print(f'final: {points_given}')

ind_list = [None] * len(points_given)
for ind in range(len(points_given)):
    curMin = 1
    for indM in range(len(points_msh)):
        if abs(points_given[ind][0]-points_msh[indM][1]) < 0.0015 and abs(points_given[ind][1]-points_msh[indM][2]) < 0.0015:
            if ind_list[ind] == None: 
                ind_list[ind] = indM
            else:
                if points_msh[indM][0] < curMin:
                    ind_list[ind] = indM
            curMin = min(curMin, points_msh[indM][0])
print("list of indices you need: ", ind_list)
#if the second index turns out to be "none", flip the angle of the rotational matrix
print(len(ind_list))