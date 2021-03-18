from abaqus import *
import testUtils
from abaqusConstants import *
import random
import math
import regionToolset
import csv
import os
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

columns = ['Cube Size','Strain','Sphere Radius','Instances','Volume Fraction','Rel Mesh Size','Elements','Nodes','Force','Lateral Disp']
with open('elasticity.csv', 'wb') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(columns)

for k in range(0,200,10):
    CUBE_SIZE = 10.0
    RADIUS = 1.0
    INSTANCES = k
    STRAIN = -0.01
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

    #Create Sphere : The sphere is created using a semicircle and a construction line to revolve it around the construction line
    sketch = model.ConstrainedSketch(name='Sketch_B', sheetSize=200.0) # for drawing the sketch ie semi circle
    sketch.ConstructionLine(point1=(0.0, -100.0), point2=(0.0, 100.0)) # construction line for revloving the semi circle
    sketch.FixedConstraint(entity=sketch.geometry[2])
    sketch.ArcByCenterEnds(center=(0.0, 0.0), direction=CLOCKWISE, point1=(0.0, -RADIUS), point2=(0.0, RADIUS))# Generating the arc
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
    model.materials['Copper'].Elastic(table=((128.0E9, 0.34), ))

    model.Material(name='Moly')
    model.materials['Moly'].Elastic(table=((329.0E9, 0.31), ))

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
        rve = assem.InstanceFromBooleanMerge(instances=[matrix]+sphereList+cutSphereList, name='rve', keepIntersections=True)
    else:
        rve = matrix
    rve_ = model.parts[rve.partName]

    # VOLUME FRACTION
    vf = 1 - matrix_.getVolume() / rve_.getVolume()

    # MESH RVE

    rve_.setMeshControls(elemShape=TET, regions=rve_.cells, technique=FREE)
    rve_.seedPart(deviationFactor=0.2, minSizeFactor=0.2, size=REL_MESH_SIZE*CUBE_SIZE)
    rve_.generateMesh()

    model.StaticStep(name='Loading', previous='Initial',
        description='Apply Load', nlgeom=OFF)

    # DEFINE SURFACES - LOAD IN Z DIRECTION

    tol = 0.01 # Extra dimensions of bounding box around face
    # Z MIN BASE
    s=rve.faces.getByBoundingBox(-tol,-tol,-tol,CUBE_SIZE+tol,CUBE_SIZE+tol,tol)
    baseSurface = assem.Surface(name='Base', side1Faces=s)
    baseSet = assem.Set(name='Base', faces=s)
    # Z MAX LOAD
    s=rve.faces.getByBoundingBox(-tol,-tol,CUBE_SIZE-tol,CUBE_SIZE+tol,CUBE_SIZE+tol,CUBE_SIZE+tol)
    loadSurface = assem.Surface(name='Load', side1Faces=s)
    loadSet = assem.Set(name='Load', faces=s)
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

    assem.regenerate()
    # BOUNDARY CONDITIONS AND COUPLING CONSTRAINTS

    model.DisplacementBC(name='Pinned Base',
        createStepName='Initial', region=baseSet, u1=UNSET, u2=UNSET, u3=SET,
        ur1=UNSET, ur2=UNSET, ur3=UNSET, amplitude=UNSET, distributionType=UNIFORM,
        fieldName='', localCsys=None)
    model.DisplacementBC(name='Pinned Side X',
        createStepName='Initial', region=sideSet1, u1=SET, u2=UNSET, u3=UNSET,
        ur1=UNSET, ur2=UNSET, ur3=UNSET, amplitude=UNSET, distributionType=UNIFORM,
        fieldName='', localCsys=None)
    model.DisplacementBC(name='Pinned Side Y',
        createStepName='Initial', region=sideSet3, u1=UNSET, u2=SET, u3=UNSET,
        ur1=UNSET, ur2=UNSET, ur3=UNSET, amplitude=UNSET, distributionType=UNIFORM,
        fieldName='', localCsys=None)
    model.DisplacementBC(name='Imposed Deformation',
        createStepName='Loading', region=loadSet, u1=UNSET, u2=UNSET, u3=CUBE_SIZE*STRAIN,
        ur1=UNSET, ur2=UNSET, ur3=UNSET, amplitude=UNSET, distributionType=UNIFORM,
        fieldName='', localCsys=None)

    v = rve.vertices.getByBoundingBox(CUBE_SIZE-tol,CUBE_SIZE-tol,CUBE_SIZE-tol,CUBE_SIZE+tol,CUBE_SIZE+tol,CUBE_SIZE+tol)
    p_max = regionToolset.Region(vertices=v)
    p_maxSet = assem.Set(vertices=v, name='p_max')
    p_maxNode = rve_.nodes.getByBoundingBox(CUBE_SIZE-tol,CUBE_SIZE-tol,CUBE_SIZE-tol,CUBE_SIZE+tol,CUBE_SIZE+tol,CUBE_SIZE+tol)[0]

    # EQUATION CONSTRAINTS - X
    i = 0
    for node in sideSet2.nodes:
        if node.label != p_maxNode.label:
            rve_.Set(name='X'+str(i), nodes=rve_.nodes[(node.label-1):(node.label)])
            model.Equation(name='X_Flat_Surface_' + str(i), terms=((1.0, rve.name+'.X' + str(i), 1), (-1.0, 'p_max', 1)))
        i += 1

    # EQUATION CONSTRAINTS - Y
    i = 0
    for node in sideSet4.nodes:
        if node.label != p_maxNode.label:
            rve_.Set(name='Y'+str(i), nodes=rve_.nodes[(node.label-1):(node.label)])
            model.Equation(name='Y_Flat_Surface_' + str(i), terms=((1.0, rve.name+'.Y' + str(i), 2), (-1.0, 'p_max', 2)))
        i += 1

    # INTEGRATE OVER SURFACES

    model.IntegratedOutputSection(name='I-Section-1',
        surface=assem.surfaces['Load'])
    model.IntegratedOutputSection(name='I-Section-2',
        surface=assem.surfaces['Base'])
    model.HistoryOutputRequest(name='H-Output-1',
        createStepName='Loading', variables=('SOF', ),
        integratedOutputSection='I-Section-1', sectionPoints=DEFAULT,
        rebar=EXCLUDE)
    model.HistoryOutputRequest(name='H-Output-2',
        createStepName='Loading', variables=('SOF', ),
        integratedOutputSection='I-Section-2', sectionPoints=DEFAULT,
        rebar=EXCLUDE)

    # SUBMIT JOB

    job = mdb.Job(name='Job-'+str(k), model=model.name, description='',
        type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0, queue=None,
        memory=90, memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True,
        explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF,
        modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='',
        scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=1,
        numGPUs=0)

    job.submit(consistencyChecking=OFF)
    job.waitForCompletion()

    # FETCH OUTPUT VARIABLES

    odb = session.openOdb(name='Job-'+str(k)+'.odb')

    HO = odb.steps['Loading'].historyRegions['Surface LOAD'].historyOutputs
    force = HO[HO.keys()[-1]].data[-1][-1]

    U = odb.steps['Loading'].frames[-1].fieldOutputs['U']
    U_p_max = U.getSubset(region=odb.rootAssembly.nodeSets['P_MAX']).values[0].data
    U_p_max = 0.5 * (U_p_max[0] + U_p_max[1]) # Average of X and Y displacement of P max

    data = [CUBE_SIZE, STRAIN, RADIUS, INSTANCES, vf, REL_MESH_SIZE, len(rve_.elements), len(rve_.nodes), force, U_p_max]

    with open('elasticity.csv', 'ab') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(data)
