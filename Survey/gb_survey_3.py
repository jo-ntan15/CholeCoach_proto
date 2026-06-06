#place the new .obj file in Student Work -> assets -> gallbladder -> lvl3
#replace just the file name in colFile and visFile
#if anything breaks for some reason dont hesistate to let me know 

import Sofa
import Sofa.Gui
import SofaRuntime
import numpy as np
import meshio
import math
import tkinter as tk
from tkinter import ttk
from multiprocessing import Process, Value, Array
import ctypes
import os


SofaRuntime.importPlugin("Sofa.Component")
SofaRuntime.importPlugin("Sofa.GL.Component")
SofaRuntime.importPlugin("CImgPlugin")

mouseExist = False
mousePosition = [0, 0, 0]
constrainedIndices = []
fixedIndices = []
constraintPlane = []
fixedPlane = [-1, 0, 0]
stiffness = Value('d',35000)
exit_flag = Value('b', False)
forceExerted = Value('d', 0.0)
testerName = Array(ctypes.c_char,100)
testerID = Array(ctypes.c_char, 100)
visual_tog = Value(ctypes.c_bool, False)
mshFile = '../assets/gallbladder/lvl3/gb_3.msh'
colFile = '../assets/gallbladder/lvl3/gb_level3.obj'
visFile = '../assets/gallbladder/lvl3/gb_level3.obj'
pngFile = '../assets/gallbladder/lvl3/frame_0188_ss.png'
points = []

taskState = [0,0,0,0,0,0]
fileName = ''
#camera var
cameraPosMax = [-0.145369, 0.0408972, 0.0802924]
cameraOrInit = [0.47652,-0.497717,-0.513223,0.521213]
cameraLookInit = [0.0619643, 0.0400911, 0.0657179]
camPosHist = [[],[]]
camOrHist = [[],[]]
camLookHist = [[],[]]
step = 0 
camRadMax = np.dot(cameraPosMax,cameraPosMax)

taskPosLock = [
    [-0.0799362,-0.029092,0.0686761],
    [-0.0864609,0.103829,0.0449736], #at around 8-9 ish minutes in
    [-0.1454,0.04089,0.08005]
]
taskOrLock = [
    [0.597186,-0.336771,-0.350834,0.637863],
    [-0.275843,0.701919,0.57237,-0.321889],
    [0.477,-0.4881,-0.5129,0.5208]
]
taskLookLock = [
    [0.0616279,0.0584224,0.0586889],
    [0.0415292,-0.000531272,0.0679091],
    [0.04109,0.04017,0.06725]
]

def intersection(lst1):
    global constrainedIndices
    return np.array([int(value) for value in lst1 if value not in constrainedIndices])
def mag(v):
    sum = 0
    for i in v:
        sum += i**2
    return math.sqrt(sum)

def createPlane(p1, p2, p3):
    p1 = np.array(p1)
    p2 = np.array(p2)
    p3 = np.array(p3)    
    v1 = p3 - p1
    v2 = p2 - p1
    normal = np.cross(v1, v2)
    D = -np.dot(normal, p1)
    return (*normal, D)

def gui(stiffness, force, testerID, visual_tog, exit_flag):
    root = tk.Tk()
    root.geometry('300x300')
    root.resizable(True, True)
    root.title('Interface')
    root.attributes('-topmost', True)
    
    # def update_stiffness(val):
    #     stiffness.value = float(val)

    def updateTester():
        cleanID = ''.join(testerIDEntry.get().split())
        testerID.value = cleanID.encode()
        if cleanID == '':
            submitLabel.config(text='ID cannot be empty')
        else: 
            submitLabel.config(text='Submitted!')

    def close():
        exit_flag.value = True
        root.destroy()
    
    def vis_switch():
        if visual_tog.value:
            togButton.config(image=off)
            visual_tog.value = False
        else:
            togButton.config(image=on)
            visual_tog.value = True
    
    testerIDFrame = ttk.Frame(root,width=250,height=40)
    testerIDFrame.pack()
    testerIDFrame.pack_propagate(False)
    testerIDLabel = ttk.Label(testerIDFrame, text='       ID: ')
    testerIDEntry = ttk.Entry(testerIDFrame)
    testerIDSub = tk.Button(testerIDFrame, command=updateTester, text = 'Submit')

    submitFrame = ttk.Frame(root)
    submitFrame.pack()
    submitLabel = ttk.Label(submitFrame,text='Please enter your ID above')
    submitLabel.pack()

    testerIDLabel.pack(side=tk.LEFT)
    testerIDEntry.pack(side=tk.LEFT, padx=5)
    testerIDSub.pack(side=tk.LEFT)

    forceFrame = ttk.Frame(root)
    forceFrame.pack(pady=20)
    force_text = ttk.Label(forceFrame, text='Force Applied = ',font=("",12))
    force_text.pack(side = tk.LEFT)
    forceLabel = ttk.Label(forceFrame, text=str(force.value), font=("",12))
    forceLabel.pack(side = tk.LEFT)
    
    on = tk.PhotoImage(file='on.png')
    off = tk.PhotoImage(file='off.png')

    togFrame = ttk.Frame(root)
    togFrame.pack(pady=10)
    togButton = ttk.Button(togFrame, image=off, command=vis_switch)
    togButton.pack()
    togBtnLab = ttk.Label(togFrame, text='Show Stress Map')
    togBtnLab.pack()

    instFrame = ttk.Frame(root)
    instFrame.pack()
    inst1 = ttk.Label(instFrame, text='Click animate and check the "real time" box.')
    inst2 = ttk.Label(instFrame, text='LMB to rotate, RMB to pan.')
    inst3 = ttk.Label(instFrame, text='Hold shift and manipulate with LMB.')
    inst1.pack()
    inst2.pack()
    inst3.pack()
    def update_force_display():
        forceLabel.config(text=f"{force.value:.2f} N")
        root.after(100, update_force_display)  # Update every 100ms

    tk.Button(root, text="Exit", command=close).pack(side=tk.BOTTOM)
    
    update_force_display()
    root.mainloop()


def visualRefresh(self):
    root = self.getContext()
    root.Gallbladder.removeChild(root.Gallbladder.getChild('Visual'))
    root.Gallbladder.addChild('Visual')
    root.Gallbladder.Visual.addObject('MeshOBJLoader', name='Surface', filename=visFile,rotation='90 0 0')
    root.Gallbladder.Visual.addObject('OglModel',name='VisualModel',src='@Surface',texturename=pngFile)
    root.Gallbladder.Visual.addObject('BarycentricMapping',name='oglMapping',input='@../dofs',output='@VisualModel')

    # root.Liver.removeChild(root.Liver.getChild('LiverVisual'))
    # root.Liver.addChild('LiverVisual')
    # root.Liver.LiverVisual.addObject('OglModel',name='VisualModel',src='@../LiverSurface',texturename='../assets/gallbladder/liver-texture-square.png')
    # root.Fat.removeObject(root.Fat.getObject('fatVisual'))
    # root.Fat.addObject('OglModel',name='fatVisual',src='@fatSurface',texturename='../assets/gallbladder/lvl3/level3_fat_updated_Material.001.png')
    root.Background.removeObject(root.Background.getObject('BgVisual'))
    root.Background.addObject('OglModel',name='BgVisual',src='@BgPlane',texturename='../assets/gallbladder/lvl3/plane_level3_updated_Material.002.png')
    [969, 315, 1253, 1220]
def activateTask(self, taskNum):
    root = self.getContext()
    global fileName
    indexLs = [739, 664,246] #highlighted index corresponding to task number
    os.makedirs(f'survey_results/{testerID.value.decode()}', exist_ok=True)
    fileName = f'survey_results/{testerID.value.decode()}/{testerID.value.decode()}_lvl3_task{taskNum+1}.txt'
    try:#[  41  144  150  151  206  207  627  628 1142 2039 2087]
        with open(fileName,'x') as out:
            out.write(f'{testerID.value.decode()} Task{taskNum+1}\n')
    except FileExistsError:
        pass
    posList = []
    for coord in points[indexLs[taskNum]]:
        posList.append(float(coord))

    root.addChild('highlight')
    root.highlight.addObject('MeshGmshLoader', name='msh',filename='../assets/gallbladder/proxySphere(120,0.005).msh',scale=2, translation=posList)
    root.highlight.addObject('TetrahedronSetTopologyContainer',name='tetraCon',src='@msh')
    root.highlight.addObject('MechanicalObject',name='HLdofs',src='@msh',template='Vec3d')
    root.highlight.addObject('UniformMass',totalMass=10000)
    root.highlight.addObject('TetrahedronFEMForceField',name='HLFEM',youngModulus='0.000001')
    root.highlight.addObject('MeshSTLLoader', name='loader', filename='../assets/gallbladder/proxySphere(120,0.005).stl',scale=2,translation=posList)
    root.highlight.addObject('OglModel', src='@loader',name='visu', texturename='../assets/gallbladder/solid_green.png')
    root.highlight.addObject('FixedProjectiveConstraint',indices=list(range(63)))
    
    
    # highlight_name = f"highlight_task{taskNum+1}"
    # root.addChild(highlight_name, name=highlight_name)

    # highlight_node = getattr(root, highlight_name)

    # highlight_node.addObject(
    #     'MeshSTLLoader',
    #     name='loader',
    #     filename='../assets/gallbladder/proxySphere(120,0.005).stl',
    #     scale=2.5,
    #     translation=posList
    # )
    # highlight_node.addObject('MeshGmshLoader', name='msh',filename='../assets/gallbladder/proxySphere(120,0.005).msh',scale=2, translation=posList)
    # highlight_node.addObject('TetrahedronSetTopologyContainer',name='tetraCon',src='@msh')
    # highlight_node.addObject(
    #     'OglModel',
    #     src='@loader',
    #     name='visu',
    #     texturename='../assets/gallbladder/solid_green.png'
    # )


    visualRefresh(self)
    Sofa.Simulation.init(root)

class Controller(Sofa.Core.Controller):
    def __init__(self, *args, **kwargs):
        Sofa.Core.Controller.__init__(self, *args, **kwargs)

    def onAnimateBeginEvent(self, event):
        root = self.getContext()  
        if visual_tog.value:
            root.dispStyle.displayFlags = "hideVisualModels showBehavior"
        else:
            root.dispStyle.displayFlags = "showVisualModels hideBehavior"
        #camera control module
        global step, camPosHist, camOrHist, camLookHist
        if not step:
            root.Camera.position.value = cameraPosMax
            root.Camera.orientation.value = cameraOrInit
            root.Camera.lookAt.value = cameraLookInit
            camPosHist = [cameraPosMax,cameraPosMax]
            camOrHist = [cameraOrInit,cameraOrInit]
            camLookHist = [cameraLookInit,cameraLookInit]
        step = 1
        camRad = np.dot(root.Camera.position.value,root.Camera.position.value)
        # print('curPos is: ',root.Camera.position.value)
        # print('curr camRad is: ',camRad)
        # print(camPosHist)
        currPos = root.Camera.position.value.tolist()
        if any(taskState):
            taskInd = taskState.index(True)
            root.Camera.position.value = taskPosLock[taskInd]
            root.Camera.orientation.value = taskOrLock[taskInd]
            root.Camera.lookAt.value = taskLookLock[taskInd]
        else:
            if camRad <= camRadMax+0.02 and currPos[0] > -0.15 and currPos[1] > -0.05 and currPos[1] < 0.15 and currPos[0] < 0.01 and currPos[2] < 0.13 and currPos[2] > 0.01:
                camPosHist.pop(0)
                camOrHist.pop(0)
                camLookHist.pop(0)
                # print('updating')
                camPosHist.append(root.Camera.position.value.tolist())
                camOrHist.append(root.Camera.orientation.value.tolist())
                camLookHist.append(root.Camera.lookAt.value.tolist())  
            else:
                root.Camera.position.value = camPosHist[0]
                root.Camera.orientation.value = camOrHist[0]
                root.Camera.lookAt.value = camLookHist[0]
        #camera end
        
        root.Gallbladder.FEM.youngModulus[0] = stiffness.value
        global mousePosition, mouseExist, forceExerted
        if not hasattr(root, "Mouse"):
            mouseExist = False
            # print('no mouse, deleting spring, ROI')
            root.Proxy.proxyObj.position = [[0,0,0]]
            if hasattr(root.Gallbladder,'spring'):
                root.Gallbladder.removeObject(root.Gallbladder.getObject('spring'))
                forceExerted.value = 0
                # print('spring deletion')
                visualRefresh(self)
                root.Gallbladder.FEM.showStressAlpha = 1
                Sofa.Simulation.init(root.Gallbladder)    

            if hasattr(root.Gallbladder,'ROI'):
                root.Gallbladder.removeObject(root.Gallbladder.getObject('ROI'))
                # print('ROI deletion')
                visualRefresh(self)
                Sofa.Simulation.init(root.Gallbladder)    

        else:
            # print('found mouse, create ROI if no ROI found AND no spring found')
            mouseExist = True
            mousePosition = root.Mouse.MousePosition.position.value[0]
            if not hasattr(root.Gallbladder,'ROI') and not hasattr(root.Gallbladder, 'spring'):
                print('ROI not found')
                root.Gallbladder.addObject('SphereROI', name='ROI', centers='@../Mouse/MousePosition.position', radii=0.005, position='@dofs.position')
                # print('made ROI and attachment proxy')
                visualRefresh(self)

            elif hasattr(root.Gallbladder,'ROI') and not hasattr(root.Gallbladder,'spring'):
                print('ROI was found, relocate for spring')
                # print('scan for nodes')
                print(root.Gallbladder.ROI.indices.value)
                indices = intersection(root.Gallbladder.ROI.indices.value)
                if len(indices) < 3: 
                    if hasattr(root.Gallbladder,'spring'):
                        root.Gallbladder.removeObject(root.Gallbladder.getObject('spring'))
                        Sofa.Simulation.init(root.Gallbladder)
                    return
                print('creating spring')
                indices2 = []
                stiffList = []
                lenlist = []
                for point in indices:
                    indices2.append(0)
                    stiffList.append(10000)   
                    lenlist.append(0.001)
        
                root.Gallbladder.addObject('StiffSpringForceField',name='spring',object2=root.Gallbladder.dofs.getLinkPath(),object1=root.Proxy.proxyObj.getLinkPath(),indices2=indices,indices1=indices2,length=lenlist,stiffness=stiffList,drawMode=0,showArrowSize=0.0002)
                if any(taskState):
                    with open(fileName, 'a') as out:
                        out.write('new spring created \n')
                visualRefresh(self)
                root.Gallbladder.FEM.showStressAlpha = 0.21
                Sofa.Simulation.init(root.Gallbladder)
                
            else:
                # print('moving spring')   
                proxyPosition = root.Proxy.proxyObj.position.value[0]
                sumIndex = 0
                forceVector = []
                forcePosition = [0,0,0]
                for i in range(3):
                    res = mousePosition[i]- proxyPosition[i]
                    forceVector.append(float(res))
                    sumIndex += res ** 2
                distance = math.sqrt(sumIndex)
                for i in range(3):
                    forceVector[i] = round(forceVector[i]/distance,2)
                    forcePosition[i] = float(round(proxyPosition[i]*1000,2))
                if distance < 0.1: 
                    forceExerted.value = round(100 * distance,2)
                if any(taskState):
                    with open(fileName, 'a') as out:
                        out.write(f"Position: {forcePosition} | Force: {forceExerted.value} | Vector: {forceVector}\n")

    def onKeypressedEvent(self, event) :
        root = self.getContext()
        global taskA, taskB, taskC, points, constraintIndices
        key = event['key']
        if key == 'L':
            root.Camera.position.value = cameraPosMax
            root.Camera.orientation.value = cameraOrInit
            root.Camera.lookAt.value = cameraLookInit 
        elif int(key) > 4:
            print('out of bounds')
        else:
            taskNum = int(key)-1
            if any(taskState):
                if taskNum != taskState.index(True):
                    print('error')
                else:
                    taskState[taskNum] = False
                    root.removeChild(root.getChild('highlight'))
                    visualRefresh(self)
                    Sofa.Simulation.init(root)                    
            else:
                taskState[taskNum] = True
                activateTask(self, taskNum)  
    
def createScene(root):
    print('creating scene')
    global constrainedIndices
    global constraintPlane
    root.dt = 0.02
    root.gravity=[0, -0.00001, 0]
    root.bbox = [[-0.04,-0.04,-0.04],[0.04,0.04,0.04]]
    print(points[739])#[739, 664,246]
    print(points[664])
    print(points[246])
    #dependencies that may or may not be necessary
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

    settings = root.addChild('Settings')
    settings.addObject('SofaDefaultPathSetting')
    app = settings.addChild('Application')
    app.addObject('BackgroundSetting', color='0 0 0',image='../assets/gallbladder/backdrop.png')
    # mouseConfig = settings.addChild('MouseConfiguration')
    # mouseConfig.addObject('AttachBodyButtonSetting', button='Left', stiffness='100')
    root.addObject('VisualStyle', name='dispStyle',displayFlags="showVisualModels hideBehavior")
    root.addObject('InteractiveCamera',name='Camera')

    root.addObject('LightManager')
    root.addObject('SpotLight', name='light1', color='1 1 1', position='-0.5 0 0', direction='1 0 0', cutoff='30', exponent='1')
    root.addObject('SpotLight', name='light2', color='1 1 1', position='0 -0.5 0', direction='0 1 0', cutoff='30', exponent='1')
    root.addObject('SpotLight', name='light3', color='1 1 1', position='0 0 -0.5', direction='0 0 1', cutoff='30', exponent='1')
    
    root.addObject('EulerImplicitSolver', rayleighStiffness='0.1', rayleighMass='0.1')
    root.addObject('CGLinearSolver', iterations='25', tolerance="1e-5", threshold="1e-5")

    root.addObject(Controller(name="proxyController"))

    gallbladder = root.addChild('Gallbladder')
    gallbladder.addObject('MeshGmshLoader',name='volMesh',filename=mshFile,rotation='0 0 0')
    gallbladder.addObject('TetrahedronSetTopologyContainer',name='volTopo',src='@volMesh')
    gallbladder.addObject('MechanicalObject',name='dofs',src='@volMesh',template="Vec3d")
    gallbladder.addObject('UniformMass',totalMass='0.04')
    gallbladder.addObject('TetrahedronFEMForceField',name='FEM',youngModulus=stiffness.value,poissonRatio='0.4',method='large',computeVonMisesStress='1',showVonMisesStressPerNodeColorMap='1')
    gallbladder.addObject('TetrahedronSetGeometryAlgorithms')
    # print(constrainedIndices)
    gallbladder.addObject('RestShapeSpringsForceField',name='conTissue',points=np.array(constrainedIndices),stiffness=12,drawSpring=1)
    gallbladder.addObject('PlaneForceField',normal='-1 0 0', d="-0.015", stiffness='1000',damping='1',showPlane='0')
    gallbladder.addObject('FixedProjectiveConstraint', indices=fixedIndices)
    # [0.98651859 0.16364923]
    # gallbladder.addObject('PlaneForceField',normal=constraintPlane[:-1],d=-1,showPlane=1)
    # gallbladder.addObject('SphereROI',name='forceROI',centers='0 0.0470 0.0560',radii=0.0075,position='@dofs.position',drawROI=1)
    # forceM = 5/len(gallbladder.forceROI.indices.value)
    # forceV = [0, forceM*0.9806, forceM*-0.1961]
    # forceV_all = []
    # for i in gallbladder.forceROI.indices.value:
    #     forceV_all.extend(forceV)
    # gallbladder.addObject('ConstantForceField', name='pull', indices=root.Gallbladder.forceROI.indices.value, forces=forceV_all)

    # global oPos
    # for point in ind_ROI:
    #     if point == None: 
    #         oPos.append(0)
    #         continue
    #     oPos.append(np.array(gallbladder.dofs.position.value[point].tolist()))
    #     posList = []
    #     for coord in points_o[point]:
    #         posList.append(float(coord))
        
    #     gallbladder.addObject('SphereROI',centers=[posList], radii=0.002,position='@dofs.position',drawROI=1)
    
    collision = gallbladder.addChild('Collision')
    collision.addObject('MeshOBJLoader',name='collisionLoader', filename=colFile, rotation='90 0 0')
    collision.addObject('TriangleSetTopologyContainer',src='@collisionLoader',name='collisionTopo')
    collision.addObject('MechanicalObject',name='collisionDOFs',src='@collisionLoader')
    collision.addObject('TriangleCollisionModel',contactStiffness='6000000',contactFriction='1')
    collision.addObject('BarycentricMapping',input="@../dofs",output='@collisionDOFs')

    visu = gallbladder.addChild('Visual')
    visu.addObject('MeshOBJLoader', name='Surface', filename=visFile,rotation='90 0 0')
    visu.addObject('OglModel', name='VisualModel', src='@Surface', texturename=pngFile)
    visu.addObject('BarycentricMapping',name='oglMapping',input='@../dofs',output='@VisualModel')
#extra stuff 
    # conn_tissue = root.addChild('ConnectiveTissue')
    # conn_tissue.addObject('MeshSTLLoader', name='LiverSurface',filename='../assets/gallbladder/connective_tissue.stl')
    # con_visu = conn_tissue.addChild('ConnTissueVisul')
    # con_visu.addObject('OglModel',name='VisualModel',src='@../LiverSurface',texturename='../assets/gallbladder/connective_tissue.png')

    proxy = root.addChild('Proxy')
    proxy.addObject('MechanicalObject',name='proxyObj',position='0 0 0')    
    proxy.addObject('UniformMass', totalMass=0.001)
    proxy.addObject('SphereCollisionModel', radius=0.007)

    # liver = root.addChild('Liver')
    # liver.addObject('MeshGmshLoader',name='LiverMesh',filename="../assets/gallbladder/liver_mesh.msh", rotation='90 0 0')
    # liver.addObject('TetrahedronSetTopologyContainer',name='liverTopo',src='@LiverMesh')
    # liver.addObject('MechanicalObject',name='liverdofs',src='@LiverMesh',template='Vec3d')
    # liver.addObject('UniformMass',totalMass=10000000)
    # meshLiver = meshio.read('../assets/gallbladder/liver_mesh.msh')
    # liverIndicesLs = []
    # liverIndicesLs.extend(range(len(meshLiver.points)))
    # liver.addObject("FixedProjectiveConstraint", indices=liverIndicesLs)
    # liver.addObject('TetrahedronFEMForceField',name='liverFEM',youngModulus='0.0000001',poissonRatio='0.4',method='large',showVonMisesStressPerNodeColorMap=1)
    # liver.addObject('MeshOBJLoader',name='LiverSurface',filename="../assets/gallbladder/liver_texturemappping_resize.obj",rotation='90 0 0', translation='0 0 0')    
    # liver_visu = liver.addChild('LiverVisual')
    # liver_visu.addObject('OglModel',name='VisualModel',src='@../LiverSurface',texturename='../assets/gallbladder/liver-texture-square.png')


    # fat = root.addChild('Fat')
    # fat.addObject('MeshOBJLoader',name='fatSurface',filename="../assets/gallbladder/lvl3/level3_fat_updated.obj",translation='0 -0.0065 0.003', rotation='-20 0 0')
    # fat.addObject('OglModel',name='fatVisual',src='@fatSurface',texturename='../assets/gallbladder/lvl3/level3_fat_updated_Material.001.png')

    #plane stuff
    background = root.addChild('Background')
    background.addObject('MeshOBJLoader', name="BgPlane", filename='../assets/gallbladder/lvl3/plane_level3_updated.obj')
    background.addObject('OglModel',name='BgVisual',src="@BgPlane",texturename='../assets/gallbladder/lvl3/plane_level3_updated_Material.002.png')
    print(points[739])

def main():
    #generating constraints
    print('running main')
    global constrainedIndices
    global constraintPlane
    global points_o, points
    mesh = meshio.read(mshFile)
    points_o = mesh.points
    surface = set()
    for cell in mesh.cells:
        if cell.type in ['triangle','quad','line']:
            surface.update(cell.data.flatten())
    
    #identifying max and min x y z
    points = mesh.points[list(surface)]
    constrainedIndices=[]
    xMinT = xMin = yMin = zMinT = zMin = 10 #xMinTotal for wall constraint, xMin for plane
    xMaxT = xMax = yMax = zMaxT = zMax = -10
    index = xMaxI = yMaxI = zMaxI = xMinI = yMinI = zMinI = 0
    for point in points:
        xMinT = min(xMinT,point[0])
        xMaxT = max(xMaxT,point[0])
        zMinT= min(zMinT,point[2])
        zMaxT = max(zMaxT,point[2])
        if point[2] < 0.1: #gating: ensuring constraint plane is created at the lower segment 
            xMin = min(xMin, point[0])
            xMax = max(xMax, point[0])
            
            yMin = min(yMin,point[1])
            yMax = max(yMax,point[1])

            zMin = min(zMin,point[2])
            zMax = max(zMax,point[2])
    # print(xMaxT,xMinT)
    xBar = xMaxT - (xMaxT - xMinT)/2
    zBar = zMaxT - (zMaxT - zMinT)/25
    # print(xBar)
    #identifying max min x y z indexes
    for point in points: 
        if point[0] >= xBar or point[2] >= zBar:
            # print(point[0],xBar)
            constrainedIndices.append(index)
        if point[0] == xMin:
            xMinI = index
            # print(index, 'xMinI update')
        if point[1] == yMin:
            yMinI = index
            # print(index, 'yMinI update')
        if point[2] == zMin: 
            zMinI = index
            # print(index, 'zMinI update')
        index += 1
    # print(points[xMinI] ,points[yMinI], points[zMinI])

    npPlane = createPlane(points[xMinI] ,points[yMinI], points[zMinI])
    constraintPlane = [float(x) for x in npPlane]
    # print(constraintPlane)

    index = 0
    for point in points:
        res = np.dot(point, constraintPlane[:-1])
        # print(res)
        if res > -0.000028: 
            if index not in constrainedIndices: 
                constrainedIndices.append(index)
                # print("plane detected and added points", index)
        
        res2 = np.dot(point, fixedPlane)
        if res2 < -0.018:
            fixedIndices.append(index)
        index += 1
    
    # print(constrainedIndices)
    # print(len(constrainedIndices))
    # print("list of fixed",fixedIndices)
    gui_process = Process(target=gui, args=(stiffness, forceExerted, testerID, visual_tog, exit_flag))
    gui_process.start()

    #starting simulation
    root = Sofa.Core.Node('root')
    createScene(root)

    Sofa.Simulation.init(root)
    
    Sofa.Gui.GUIManager.Init("myscene", "qt")
    Sofa.Gui.GUIManager.createGUI(root, __file__)

    Sofa.Gui.GUIManager.SetDimension(2000, 1500)
    
    print("scene graph loaded")
    Sofa.Gui.GUIManager.MainLoop(root)
    Sofa.Gui.GUIManager.closeGUI()

if __name__ == '__main__':
    main()