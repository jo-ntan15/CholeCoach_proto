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


SofaRuntime.importPlugin("Sofa.Component")
SofaRuntime.importPlugin("Sofa.GL.Component")
SofaRuntime.importPlugin("CImgPlugin")

mouseExist = False
mousePosition = [0, 0, 0]
constrainedIndices = []
fixedIndices = []
constraintPlane = []
fixedPlane = [-1, 0, 0]

stiffness = Value('d', 60000)
exit_flag = Value('b', False)
forceExerted = Value('d', 0.0)

testerName = Array(ctypes.c_char,100)
testerID = Array(ctypes.c_char, 100)
visual_tog = Value(ctypes.c_bool, False)

mshFile = '../assets/gallbladder/combined_ct_texture_2_volume.msh'
colFile = '../assets/gallbladder/reduced.obj'
visFile = '../assets/gallbladder/demo_texturemap.obj'
pngFile = '../assets/gallbladder/gb_texture2.png'
points = []

taskState = False
fileName = ''
#camera var
cameraPosMax = [-0.1464, 0.01742, 0.05086]
cameraOrInit = [0.443,-0.542,-0.4544,0.549]
cameraLookInit = [2.86338, 0.00071726, -0.019076]
camPosHist = [[],[]]
camOrHist = [[],[]]
camLookHist = [[],[]]
step = 0 
camRadMax = np.dot(cameraPosMax,cameraPosMax)

def intersection(lst1):
    global constrainedIndices
    return np.array([int(value) for value in lst1 if value not in constrainedIndices])

def createPlane(p1, p2, p3):
    p1 = np.array(p1)
    p2 = np.array(p2)
    p3 = np.array(p3)    
    v1 = p3 - p1
    v2 = p2 - p1
    normal = np.cross(v1, v2)
    D = -np.dot(normal, p1)
    return (*normal, D)

def gui(stiffness, force, testerName, testerID,visual_tog, exit_flag):
    root = tk.Tk()
    root.geometry('300x300')
    root.resizable(True, True)
    root.title('Interface')
    root.attributes('-topmost', True)
    
    def update_stiffness(val):
        stiffness.value = float(val)

    # def updateTester():
    #     testerName.value = testerNameEntry.get().encode()
    #     testerID.value = testerIDEntry.get().encode()

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

    sliderFrame = ttk.Frame(root)
    sliderFrame.pack(pady=10)

    slider = tk.Scale(
        sliderFrame, from_=10000, to=100000, resolution=10000,
        orient="horizontal",command=update_stiffness, length=200
    )
    slider.set(60000)
    start_label = ttk.Label(sliderFrame, text="Low")
    start_label.pack(side=tk.LEFT, padx=5)

    slider.pack(side=tk.LEFT, fill="x")
    end_label = ttk.Label(sliderFrame, text="High")
    end_label.pack(side=tk.LEFT, padx=5)
    
    # testerNameFrame = ttk.Frame(root,width=250,height=40)
    # testerNameFrame.pack()
    # testerNameFrame.pack_propagate(False)
    # testerNameLabel = ttk.Label(testerNameFrame, text='Name: ')
    # testerNameEntry = ttk.Entry(testerNameFrame)
    
    # testerNameLabel.pack(side=tk.LEFT)
    # testerNameEntry.pack(side=tk.LEFT, padx=5)
    
    # testerIDFrame = ttk.Frame(root,width=250,height=40)
    # testerIDFrame.pack()
    # testerIDFrame.pack_propagate(False)
    # testerIDLabel = ttk.Label(testerIDFrame, text='       ID: ')
    # testerIDEntry = ttk.Entry(testerIDFrame)
    # testerIDSub = tk.Button(testerIDFrame, command=updateTester, text = 'Submit')

    # testerIDLabel.pack(side=tk.LEFT)
    # testerIDEntry.pack(side=tk.LEFT, padx=5)
    # testerIDSub.pack(side=tk.LEFT)

    forceFrame = ttk.Frame(root)
    forceFrame.pack(pady=20)
    force_text = ttk.Label(forceFrame, text='Force Applied = ',font=("Segoe UI",12))
    force_text.pack(side = tk.LEFT)
    forceLabel = ttk.Label(forceFrame, text=str(force.value), font=("Segoe UI",12))
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
    inst4 = ttk.Label(instFrame, text='Change model stiffness with slider.')
    inst1.pack()
    inst2.pack()
    inst3.pack()
    inst4.pack()
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
    root.Gallbladder.Visual.addObject('OglModel',name='VisualModel',src='@Surface',texturename=pngFile,color=[1.0, 1.0, 1.0, 1.0],material='textured 1 0 0 0 1')
    root.Gallbladder.Visual.addObject('BarycentricMapping',name='oglMapping',input='@../dofs',output='@VisualModel')
    # root.Liver.removeChild(root.Liver.getChild('LiverVisual'))
    # root.Liver.addChild('LiverVisual')
    # root.Liver.LiverVisual.addObject('OglModel',name='VisualModel',src='@../LiverSurface',texturename='../assets/gallbladder/liver-texture-square.png')
    # Sofa.Simulation.init(root.Liver)
    root.Background.removeObject(root.Background.getObject('BgVisual'))
    root.Background.addObject('OglModel',name='BgVisual',src='@BgPlane',texturename='../assets/gallbladder/plane_demo_updated_Material.002.png')

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
        global step, camPosHist, camOrHist, camLookHist, taskState
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
        if taskState:
            pass
            # taskInd = taskState.index(True)
            # root.Camera.position.value = taskPosLock[taskInd]
            # root.Camera.orientation.value = taskOrLock[taskInd]
            # root.Camera.lookAt.value = taskLookLock[taskInd]
        else:
            if camRad <= camRadMax+0.01 and abs(currPos[1]) < 0.25 and currPos[0] < -0.02:
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
        # print('cam done')
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
            # print('no mouse end')
        else:
            print('found mouse, create ROI if no ROI found AND no spring found')
            mouseExist = True
            mousePosition = root.Mouse.MousePosition.position.value[0]
            if not hasattr(root.Gallbladder,'ROI') and not hasattr(root.Gallbladder, 'spring'):
                # print('ROI not found')
                root.Gallbladder.addObject('SphereROI', name='ROI', centers='@../Mouse/MousePosition.position', radii=0.005, position='@dofs.position')
                # print('made ROI and attachment proxy')
                visualRefresh(self)

            elif hasattr(root.Gallbladder,'ROI') and not hasattr(root.Gallbladder,'spring'):
                # print('ROI was found, relocate for spring')
                # print('scan for nodes')
                print(root.Gallbladder.ROI.indices.value)
                indices = intersection(root.Gallbladder.ROI.indices.value)
                if len(indices) < 5: 
                    if hasattr(root.Gallbladder,'spring'):
                        root.Gallbladder.removeObject(root.Gallbladder.getObject('spring'))
                        Sofa.Simulation.init(root.Gallbladder)
                    return
                # print('creating spring')
                indices2 = []
                stiffList = []
                lenlist = []
                for point in indices:
                    indices2.append(0)
                    stiffList.append(10000)   
                    lenlist.append(0.001)
        
                root.Gallbladder.addObject('StiffSpringForceField',name='spring',object2=root.Gallbladder.dofs.getLinkPath(),object1=root.Proxy.proxyObj.getLinkPath(),indices2=indices,indices1=indices2,length=lenlist,stiffness=stiffList,drawMode=0,showArrowSize=0.0002)
                
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
                # if any(taskState):
                #     with open(fileName, 'a') as out:
                #         out.write(f"Position: {forcePosition} | Force: {forceExerted.value} | Vector: {forceVector}\n")

        # print('end')
    def onKeypressedEvent(self, event) :
        root = self.getContext()
        global taskA, taskB, taskC, points, constraintIndices, taskState
        key = event['key']
        if key == 'L':
            root.Camera.position.value = cameraPosMax
            root.Camera.orientation.value = cameraOrInit
            root.Camera.lookAt.value = cameraLookInit 
        elif key == '1':
            if taskState:
                # print('end task')
                taskState = False
                root.removeChild(root.getChild('highlight'))
                visualRefresh(self)
                Sofa.Simulation.init(root)
            else:
                taskState = True
                posList = []
                for coord in points[69]:
                    posList.append(float(coord))

                root.addChild('highlight')
                root.highlight.addObject('MeshGmshLoader', name='msh',filename='../assets/gallbladder/proxySphere(120,0.005).msh',scale=2, translation=posList)
                root.highlight.addObject('TetrahedronSetTopologyContainer',name='tetraCon',src='@msh')
                root.highlight.addObject('MechanicalObject',name='HLdofs',src='@msh',template='Vec3d')
                root.highlight.addObject('UniformMass',totalMass=10000)
                root.highlight.addObject('TetrahedronFEMForceField',name='HLFEM',youngModulus='0.000001')
                # root.highlight.addObject('MeshSTLLoader', name='loader', filename='../assets/gallbladder/proxySphere(120,0.005).stl',scale=2,translation=posList)
                root.highlight.addObject('OglModel', src='@msh',name='visu', texturename='../assets/gallbladder/solid_green.png')
                root.highlight.addObject('FixedProjectiveConstraint',indices=list(range(63)))
                visualRefresh(self)
                Sofa.Simulation.init(root)
                
        # print('keypress fin')
        # elif int(key) > 6:
        #     print('out of bounds')
        # else:
        #     taskNum = int(key)-1
        #     if any(taskState):
        #         if taskNum != taskState.index(True):
        #             print('error')
        #         else:
        #             taskState[taskNum] = False
        #             root.removeChild(root.getChild('highlight'))
        #             visualRefresh(self)
        #             Sofa.Simulation.init(root)                    
        #     else:
        #         taskState[taskNum] = True
        #         activateTask(self, taskNum)  
                
    
def createScene(root):
    print('creating scene')
    global constrainedIndices
    global constraintPlane
    root.dt = 0.02
    root.gravity=[0, -0.00001, 0]
    root.bbox = [[-0.04,-0.04,-0.04],[0.04,0.04,0.04]]

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
    root.addObject('VisualStyle',name='dispStyle',displayFlags='showVisualModels hideBehavior')
    settings = root.addChild('Settings')
    settings.addObject('SofaDefaultPathSetting')
    app = settings.addChild('Application')
    app.addObject('BackgroundSetting', color='0 0 0',image='../assets/gallbladder/backdrop.png')
    # mouseConfig = settings.addChild('MouseConfiguration')
    # mouseConfig.addObject('AttachBodyButtonSetting', button='Left', stiffness='100')

    root.addObject('LightManager')
    root.addObject('SpotLight', name='light1', color='1 1 1', position='-0.5 0 0', direction='1 0 0', cutoff='30', exponent='1')
    root.addObject('SpotLight', name='light2', color='1 1 1', position='0 -0.5 0', direction='0 1 0', cutoff='30', exponent='1')
    root.addObject('SpotLight', name='light3', color='1 1 1', position='0 0 -0.5', direction='0 0 1', cutoff='30', exponent='1')
    
    root.addObject('InteractiveCamera',name='Camera')

    root.addObject('EulerImplicitSolver', rayleighStiffness='0.1', rayleighMass='0.1')
    root.addObject('CGLinearSolver', iterations='25', tolerance="1e-5", threshold="1e-5")

    root.addObject(Controller(name="proxyController"))

    gallbladder = root.addChild('Gallbladder')
    gallbladder.addObject('MeshGmshLoader',name='volMesh',filename=mshFile,rotation='0 0 0')
    gallbladder.addObject('TetrahedronSetTopologyContainer',name='volTopo',src='@volMesh')
    gallbladder.addObject('MechanicalObject',name='dofs',src='@volMesh',template="Vec3d")
    gallbladder.addObject('UniformMass',totalMass='0.04')
    gallbladder.addObject('TetrahedronFEMForceField',name='FEM',youngModulus='60000',poissonRatio='0.4',method='large',computeVonMisesStress='1',showVonMisesStressPerNodeColorMap='1')
    gallbladder.addObject('TetrahedronSetGeometryAlgorithms')
    print(constrainedIndices)
    gallbladder.addObject('RestShapeSpringsForceField',name='conTissue',points=np.array(constrainedIndices),stiffness=44,drawSpring=1)
    gallbladder.addObject('PlaneForceField',normal='-1 0 0', d="-0.015", stiffness='1000',damping='1',showPlane='0')
    gallbladder.addObject('FixedProjectiveConstraint', indices=fixedIndices)
    collision = gallbladder.addChild('Collision')
    collision.addObject('MeshOBJLoader',name='collisionLoader', filename=colFile)
    collision.addObject('TriangleSetTopologyContainer',src='@collisionLoader',name='collisionTopo')
    collision.addObject('MechanicalObject',name='collisionDOFs',src='@collisionLoader')
    collision.addObject('TriangleCollisionModel',contactStiffness='6000000',contactFriction='1')
    collision.addObject('BarycentricMapping',input="@../dofs",output='@collisionDOFs')

    visu = gallbladder.addChild('Visual')
    visu.addObject('MeshOBJLoader', name='Surface', filename=visFile,rotation='90 0 0',)
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
    proxy.addObject('SphereCollisionModel', radius=0.005)

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
    background = root.addChild('Background')
    background.addObject('MeshOBJLoader', name="BgPlane", filename='../assets/gallbladder/plane_demo_updated.obj',rotation='85 0 0', translation='-0.01 0 0')
    background.addObject('OglModel',name='BgVisual',src="@BgPlane",texturename='../assets/gallbladder/plane_demo_updated_Material.002.png')

def main():
    #generating constraints
    print('running main')
    global constrainedIndices
    global constraintPlane
    global points
    mesh = meshio.read(mshFile)
    points = mesh.points
    surface = set()
    for cell in mesh.cells:
        if cell.type in ['triangle','quad','line']:
            surface.update(cell.data.flatten())
    
    #identifying max and min x y z
    points = mesh.points[list(surface)]
    constrainedIndices=[]
    xMin = yMin = zMin = 10
    xMax = yMax = zMax = -10
    index = xMaxI = yMaxI = zMaxI = xMinI = yMinI = zMinI = 0
    for point in points:
        xMin = min(xMin,point[0])
        xMax = max(xMax,point[0])

        yMin = min(yMin,point[1])
        yMax = max(yMax,point[1])

        zMin = min(zMin,point[2])
        zMax = max(zMax,point[2])
    xBar = xMax - (xMax - xMin)/2.4
    print(xBar)
    #identifying max min x y z indexes
    for point in points: 
        if point[0] >= xBar:
            constrainedIndices.append(index)
        if point[0] == xMin and point[2] < 0.05: 
            xMinI = index
            print(index, 'xMinI update')
        elif point[1] == yMin and point[2] < 0.05:
            yMinI = index
            print(index, 'yMinI update')
        elif point[2] == zMin: 
            zMinI = index
            print(index, 'zMinI update')
        index += 1
    print(points[xMinI] ,points[yMinI], points[zMinI])

    npPlane = createPlane(points[xMinI] ,points[yMinI], points[zMinI])
    constraintPlane = [float(x) for x in npPlane]

    index = 0
    for point in points:
        res = np.dot(point, constraintPlane[:-1])
        if res > 0.00001: 
            if index not in constrainedIndices: 
                constrainedIndices.append(index)
                print(index)
        res2 = np.dot(point, fixedPlane)
        if res2 < -0.015:
            fixedIndices.append(index)
        index += 1
    
    print(constrainedIndices)
    print(len(constrainedIndices))
    print("list of fixed",fixedIndices)
    gui_process = Process(target=gui, args=(stiffness, forceExerted, testerName, testerID, visual_tog,exit_flag))
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