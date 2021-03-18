from abaqus import *
import testUtils
from abaqusConstants import *
import random
import math
import csv
import os
import mesh
testUtils.setBackwardCompatibility()


def dist(a, b):
    return math.sqrt((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2)


def position():
    successful = False
    tries = 0
    while not successful:
        # Trial location (sphere entirely within cube)
        guess_position = (CUBE_SIZE * random.random(), CUBE_SIZE * random.random(), CUBE_SIZE * random.random())
        successful = True
        # Check sphere will not intersect with existing spheres
        for s in sphereList:
            if dist(s.getTranslation(), guess_position) <= 2 * RADIUS:
                successful = False

                break
    return guess_position


# CLEAR TEMPORARY FILES
path = 'C:/temp'
for f in os.listdir(path):
    try:
        os.remove(os.path.join(path, f))
    except:
        pass

columns = ['Cube Size','Temp Diff','Sphere Radius','Instances','Volume Fraction','Rel Mesh Size','Elements','Nodes','Heat Flux']
with open('thermal.csv', 'wb') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(columns)

for k in range(0,140,10):
    CUBE_SIZE = 10.0
    RADIUS = 1.0
    INSTANCES = k
    TEMP_DIFF = 100
    REL_MESH_SIZE = 0.05

    model = mdb.Model(name='Particle_Composite_'+str(k))
    assem = model.rootAssembly

    # CREATE CUBE

    sketch = model.Sketch(name='Sketch_A', sheetSize=200.0)
    sketch.rectangle(point1=(0.0, 0.0), point2=(CUBE_SIZE, CUBE_SIZE))

    cube_ = model.Part(dimensionality=THREE_D, name='Cube', type=DEFORMABLE_BODY) # making it 3d
    cube_.BaseSolidExtrude(sketch=sketch, depth=CUBE_SIZE)
    del sketch

    cube = assem.Instance(dependent=ON, name='Cube', part=cube_)

    # Create Big Cube (for cutting spheres later)
    sketch = model.Sketch(name='Sketch_A', sheetSize=200.0)
    sketch.rectangle(point1=(-RADIUS, -RADIUS), point2=(CUBE_SIZE+RADIUS, CUBE_SIZE+RADIUS))

    big_cube_ = model.Part(dimensionality=THREE_D, name='Big Cube', type=DEFORMABLE_BODY) # making it 3d
    big_cube_.BaseSolidExtrude(sketch=sketch, depth=CUBE_SIZE + 2*RADIUS)
    del sketch

    big_cube = assem.Instance(dependent=ON, name='Big Cube', part=big_cube_)
    assem.translate(instanceList=('Big Cube',),vector=(0.0, 0.0, -RADIUS))

    outer_shell = assem.InstanceFromBooleanCut(cuttingInstances=(assem.instances['Cube'],), instanceToBeCut=big_cube, name='Outer Shell')

    assem.resumeAllFeatures()

    # Create Sphere : The sphere is created using a semicircle and a construction line to revolve it around the construction line
    sketch = model.ConstrainedSketch(name='Sketch_B', sheetSize=200.0) # for drawing the sketch ie semi circle
    sketch.ConstructionLine(point1=(0.0, -100.0), point2=(0.0, 100.0)) # construction line for revolving the semicircle
    sketch.FixedConstraint(entity=sketch.geometry[2])
    sketch.ArcByCenterEnds(center=(0.0, 0.0), direction=CLOCKWISE, point1=(0.0, -RADIUS), point2=(0.0, RADIUS))
    sketch.Line(point1=(0.0, -RADIUS), point2=(0.0, RADIUS)) #Closing the arc to create semi circle
    sketch.VerticalConstraint(entity=sketch.geometry[4])
    sphere = model.Part(dimensionality=THREE_D, name='Sphere', type=DEFORMABLE_BODY)# Setting the type and dimensionality of the sketch ie it should be 3d
    sphere.BaseSolidRevolve(angle=360.0, flipRevolveDirection=OFF, sketch=sketch) # Revolving the geometry
    del sketch

    sphereList = []

    # Instantiate Sphere
    for i in range(INSTANCES):

        InstanceName = 'Sphere ' + str(i)
        s = assem.Instance(dependent=ON, name=InstanceName, part=sphere)
        # Translate Instance of Sphere: the sphere instance is translated to the position generated above
        assem.translate(instanceList=(InstanceName,),vector=position())  # translate a copy of the sphere created initially
        sphereList.append(s)

    # Matrix is produced by cutting away spheres from the cube
    if sphereList:
        matrix = assem.InstanceFromBooleanCut(cuttingInstances=(sphereList), instanceToBeCut=cube, name='Matrix')
        assem.resumeAllFeatures() # Spheres were suppressed after cut
        assem.deleteFeatures(('Cube', 'Big Cube'))  # Remove temporary cube
        del model.parts[cube_.name]
    else:
        matrix = cube
        assem.deleteFeatures(('Big Cube',))

    matrix_ = matrix.part
    del model.parts[big_cube_.name]

    # CUT SPHERES

    cutSphereList = []
    temp = [] # Store spheres to be that have been cut temporarily
    for s in sphereList:
        try:
            c = assem.InstanceFromBooleanCut(cuttingInstances=(outer_shell, ), instanceToBeCut=s, name='Cut '+s.name)
            assem.resumeAllFeatures()
            cutSphereList.append(c)
            assem.deleteFeatures((s.name,)) # Delete protruding sphere
            temp.append(s)

        except AbaqusException as message:
            pass

    del assem.instances[outer_shell.name]
    for sphere_to_remove in temp:
        sphereList.remove(sphere_to_remove)

    # ASSIGN MATERIALS

    model.Material(name='Copper')
    model.materials['Copper'].Conductivity(table=((401.0, ), ))

    model.Material(name='Moly')
    model.materials['Moly'].Conductivity(table=((138.0, ), ))

    model.HomogeneousSolidSection(material='Copper', name='Matrix_Section', thickness=None)
    model.HomogeneousSolidSection(material='Moly', name='Particles_Section', thickness=None)

    matrix_.Set(cells=matrix_.cells.getSequenceFromMask(('[#1 ]',), ), name='Matrix_Set')
    matrix_.SectionAssignment(offset=0.0,
        offsetField='', offsetType=MIDDLE_SURFACE, region=
        matrix_.sets['Matrix_Set'], sectionName=
        'Matrix_Section', thicknessAssignment=FROM_SECTION)

    sphere.Set(cells=sphere.cells.getSequenceFromMask(('[#1 ]',), ), name='Sphere_Set')
    sphere.SectionAssignment(offset=0.0,
        offsetField='', offsetType=MIDDLE_SURFACE, region=
        sphere.sets['Sphere_Set'], sectionName=
        'Particles_Section', thicknessAssignment=FROM_SECTION)

    for cut_sphere in cutSphereList:
        c = cut_sphere.part
        t =c.Set(cells=c.cells.getSequenceFromMask(('[#1 ]',), ), name='Cut_Sphere_Set')
        c.SectionAssignment(offset=0.0,
                                 offsetField='', offsetType=MIDDLE_SURFACE, region=
                                 t, sectionName=
                                 'Particles_Section', thicknessAssignment=FROM_SECTION)

    # CREATE REPRESENTATIVE VOLUME ELEMENT

    if sphereList or cutSphereList:
        rve = assem.InstanceFromBooleanMerge(instances=[matrix] + sphereList + cutSphereList, name='rve',
                                             keepIntersections=True)
    else:
        rve = matrix
    rve_ = model.parts[rve.partName]

    # VOLUME FRACTION

    vf = 1 - matrix_.getVolume() / rve_.getVolume()

    # MESH RVE
    elemType = mesh.ElemType(elemCode=DC3D10, elemLibrary=STANDARD)
    rve_.setElementType(regions=(rve_.cells,), elemTypes=(elemType,))
    rve_.setMeshControls(elemShape=TET, regions=rve_.cells, technique=FREE)
    rve_.seedPart(deviationFactor=0.2, minSizeFactor=0.2, size=0.4)
    rve_.generateMesh()

    model.HeatTransferStep(name='Heating',
        previous='Initial', response=STEADY_STATE, amplitude=RAMP)

    # DEFINE SURFACES - HEATING ON Z MAX FACE

    tol = 0.01 # Extra dimensions of bounding box around face
    # Z MIN BASE
    s=rve.faces.getByBoundingBox(-tol,-tol,-tol,CUBE_SIZE+tol,CUBE_SIZE+tol,tol)
    baseSurface = assem.Surface(name='Base', side1Faces=s)
    baseSet = assem.Set(name='Base', faces=s)
    # Z MAX HEATING
    s=rve.faces.getByBoundingBox(-tol,-tol,CUBE_SIZE-tol,CUBE_SIZE+tol,CUBE_SIZE+tol,CUBE_SIZE+tol)
    heatedSurface = assem.Surface(name='Heated', side1Faces=s)
    heatedSet = assem.Set(name='Heated', faces=s)
    # X MIN
    s=rve.faces.getByBoundingBox(-tol,-tol,-tol,tol,CUBE_SIZE+tol,CUBE_SIZE+tol)
    sideSurface1 = assem.Surface(name='Side1', side1Faces=s)
    sideSet1 = assem.Set(name='Side1', faces=s)
    # X MAX
    s=rve.faces.getByBoundingBox(CUBE_SIZE-tol,-tol,-tol,CUBE_SIZE+tol,CUBE_SIZE+tol,CUBE_SIZE+tol)
    sideSurface2 = assem.Surface(name='Side2', side1Faces=s)
    sideSet2 = assem.Set(name='Side2', faces=s)
    # Y MIN
    s=rve.faces.getByBoundingBox(-tol,-tol,-tol,CUBE_SIZE+tol,tol,CUBE_SIZE+tol)
    sideSurface3 = assem.Surface(name='Side3', side1Faces=s)
    sideSet3 = assem.Set(name='Side3', faces=s)
    # Y MAX
    s=rve.faces.getByBoundingBox(-tol,CUBE_SIZE-tol,-tol,CUBE_SIZE+tol,CUBE_SIZE+tol,CUBE_SIZE+tol)
    sideSurface4 = assem.Surface(name='Side4', side1Faces=s)
    sideSet4 = assem.Set(name='Side4', faces=s)

    # BOUNDARY CONDITIONS AND COUPLING CONSTRAINTS

    model.TemperatureBC(name='Fixed temp base',
        createStepName='Initial', region=baseSet, distributionType=UNIFORM,
        fieldName='', magnitude=0.0)
    model.TemperatureBC(name='Apply heat',
        createStepName='Heating', region=heatedSet, fixed=OFF,
        distributionType=UNIFORM, fieldName='', magnitude=TEMP_DIFF, amplitude=UNSET)

    # INTEGRATE OVER SURFACES

    model.IntegratedOutputSection(name='I-Section-1',
        surface=assem.surfaces['Heated'])
    model.IntegratedOutputSection(name='I-Section-2',
        surface=assem.surfaces['Base'])
    model.HistoryOutputRequest(name='H-Output-1',
        createStepName='Heating', variables=('SOH', ),
        integratedOutputSection='I-Section-1', sectionPoints=DEFAULT,
        rebar=EXCLUDE)
    model.HistoryOutputRequest(name='H-Output-2',
        createStepName='Heating', variables=('SOH', ),
        integratedOutputSection='I-Section-2', sectionPoints=DEFAULT,
        rebar=EXCLUDE)

    # SUBMIT JOB

    job = mdb.Job(name='Job-' + str(k), model=model.name, description='',
                  type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0, queue=None,
                  memory=90, memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True,
                  explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF,
                  modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='',
                  scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=1,
                  numGPUs=0)

    job.submit(consistencyChecking=OFF)
    job.waitForCompletion()

    odb = session.openOdb(name='Job-'+str(k)+'.odb')

    HO = odb.steps['Heating'].historyRegions['Surface HEATED'].historyOutputs
    heatFlux = HO[HO.keys()[-1]].data[-1][-1]

    data = [CUBE_SIZE, TEMP_DIFF, RADIUS, INSTANCES, vf, REL_MESH_SIZE, len(rve_.elements), len(rve_.nodes), heatFlux]

    with open('thermal.csv', 'ab') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(data)