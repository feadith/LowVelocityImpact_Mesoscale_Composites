# library of npython
import os

# Save by pulungd on 2018_10_15-09.17.57; build 6.14-5 2015_08_18-17.37.49 135153
from abaqus import *
from abaqusConstants import *
from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *
#
import __main__

def createImpactorBall(BallRad,partName):
    import section
    import regionToolset
    import displayGroupMdbToolset as dgm
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import optimization
    import job
    import sketch
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    import connectorBehavior
    s1 = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
        sheetSize=200.0)
    g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
    s1.setPrimaryObject(option=STANDALONE)
    s1.ConstructionLine(point1=(0.0, -100.0), point2=(0.0, 100.0))
    s1.FixedConstraint(entity=g[2])
    s1.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(0.0, -BallRad))
    s1.CoincidentConstraint(entity1=v[1], entity2=g[2], addUndoState=False)
    s1.CoincidentConstraint(entity1=v[0], entity2=g[2], addUndoState=False)
    s1.Line(point1=(0.0, -BallRad), point2=(0.0, BallRad))
    s1.VerticalConstraint(entity=g[4], addUndoState=False)
    s1.PerpendicularConstraint(entity1=g[3], entity2=g[4], addUndoState=False)
    s1.CoincidentConstraint(entity1=v[2], entity2=g[2], addUndoState=False)
    s1.autoTrimCurve(curve1=g[3], point1=(-BallRad, -0.0))
    p = mdb.models['Model-1'].Part(name=partName, dimensionality=THREE_D, 
        type=DISCRETE_RIGID_SURFACE)
    p = mdb.models['Model-1'].parts[partName]
    p.BaseShellRevolve(sketch=s1, angle=360.0, flipRevolveDirection=OFF)
    s1.unsetPrimaryObject()
    session.viewports['Viewport: 1'].setValues(displayedObject=p)
    del mdb.models['Model-1'].sketches['__profile__']
    v2, e1, d2, n1 = p.vertices, p.edges, p.datums, p.nodes
    p.ReferencePoint(point=v2[0])
    r = p.referencePoints
    refPoints=(r[2], )
    p.Set(referencePoints=refPoints, name='Set-RP-'+partName)


def createPinHolder(pinRad,pinHeight,partName):
    import section
    import regionToolset
    import displayGroupMdbToolset as dgm
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import optimization
    import job
    import sketch
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    import connectorBehavior
    s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
        sheetSize=200.0)
    g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=STANDALONE)
    s.CircleByCenterPerimeter(center=(0.0, 0.0), point1=(0.0, -pinRad))
    p = mdb.models['Model-1'].Part(name=partName, dimensionality=THREE_D, 
        type=DISCRETE_RIGID_SURFACE)
    p = mdb.models['Model-1'].parts[partName]
    p.BaseSolidExtrude(sketch=s, depth=pinHeight)
    s.unsetPrimaryObject()
    session.viewports['Viewport: 1'].setValues(displayedObject=p)
    del mdb.models['Model-1'].sketches['__profile__']
    c1 = p.cells
    p.RemoveCells(cellList = c1[0:1])
    v1, e, d1, n = p.vertices, p.edges, p.datums, p.nodes
    p.ReferencePoint(point=p.InterestingPoint(edge=e[0], rule=CENTER))
#

def createBaseFixture(lx,ly,cutoutx,cutouty,partName):
    import section
    import regionToolset
    import displayGroupMdbToolset as dgm
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import optimization
    import job
    import sketch
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    import connectorBehavior
    s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', 
        sheetSize=200.0)
    g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
    s.setPrimaryObject(option=STANDALONE)
    s.rectangle(point1=(0.0, 0.0), point2=(lx, ly))
    s.rectangle(point1=((lx-cutoutx)/2.0, (ly-cutouty)/2.0), point2=((lx + cutoutx)/2.0, (ly + cutouty)/2.0))
    p = mdb.models['Model-1'].Part(name=partName, 
        dimensionality=THREE_D, type=DISCRETE_RIGID_SURFACE)
    p = mdb.models['Model-1'].parts[partName]
    p.BaseShell(sketch=s)
    s.unsetPrimaryObject()
    p.ReferencePoint(point= p.vertices.findAt((0.0, 0.0, 0.0), ))
    session.viewports['Viewport: 1'].setValues(displayedObject=p)
    del mdb.models['Model-1'].sketches['__profile__']
#

def MakeDatumbyOffset(myoffset,xref,yref,zref):
	import section
	import regionToolset
	import displayGroupMdbToolset as dgm
	import part
	import material
	import assembly
	import step
	import interaction
	import load
	import mesh
	import optimization
	import job
	import sketch
	import visualization
	import xyPlot
	import displayGroupOdbToolset as dgo
	import connectorBehavior
	p = mdb.models['Model-1'].parts['Part-1']
	f = p.faces
	p.DatumPlaneByOffset(plane=f.findAt(coordinates=(xref,yref, 
		zref)), flip=SIDE2, offset=myoffset)

def AssignStackDir(xloc,yloc,zloc,xref,yref,zref):
    import section
    import regionToolset
    import displayGroupMdbToolset as dgm
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import optimization
    import job
    import sketch
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    import connectorBehavior
    p = mdb.models['Model-1'].parts['Part-1']
    c = p.cells
    pickedCells = c.findAt(((xloc,yloc,zloc), ))
    f1 = p.faces
    p.assignStackDirection(referenceRegion=f1.findAt(coordinates=(xref, 
        yref, zref)), cells=pickedCells)

def AssignElType(xloc,yloc,zloc):
    import section
    import regionToolset
    import displayGroupMdbToolset as dgm
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import optimization
    import job
    import sketch
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    import connectorBehavior
    p = mdb.models['Model-1'].parts['Part-1']
    elemType1 = mesh.ElemType(elemCode=COH3D8, elemLibrary=STANDARD, 
        elemDeletion=ON, viscosity=0.0001, maxDegradation=0.999)
    elemType2 = mesh.ElemType(elemCode=COH3D6, elemLibrary=STANDARD)
    elemType3 = mesh.ElemType(elemCode=UNKNOWN_TET, elemLibrary=STANDARD)
    c = p.cells
    cells = c.findAt(((xloc, yloc, zloc), ))
    pickedRegions =(cells, )
    p.setElementType(regions=pickedRegions, elemTypes=(elemType1, elemType2, 
        elemType3))


def AssignMatl(xloc,yloc,zloc,mySetName,mySection):
	import section
	import regionToolset
	import displayGroupMdbToolset as dgm
	import part
	import material
	import assembly
	import step
	import interaction
	import load
	import mesh
	import optimization
	import job
	import sketch
	import visualization
	import xyPlot
	import displayGroupOdbToolset as dgo
	import connectorBehavior
	p = mdb.models['Model-1'].parts['Part-1']
	c = p.cells
	cells = c.findAt(((xloc, yloc, zloc), ))
	region = p.Set(cells=cells, name=mySetName)
	p.SectionAssignment(region=region, sectionName=mySection, offset=0.0, 
		offsetType=MIDDLE_SURFACE, offsetField='', 
		thicknessAssignment=FROM_SECTION)


def PartitionIt(zheight,PartitionNo):
	import section
	import regionToolset
	import displayGroupMdbToolset as dgm
	import part
	import material
	import assembly
	import step
	import interaction
	import load
	import mesh
	import optimization
	import job
	import sketch
	import visualization
	import xyPlot
	import displayGroupOdbToolset as dgo
	import connectorBehavior
	p = mdb.models['Model-1'].parts['Part-1']
	c = p.cells
	pickedCells = c.findAt(((10.0, 10.0, zheight), ))
	d = p.datums
	p.PartitionCellByDatumPlane(datumPlane=d[PartitionNo], cells=pickedCells)

def AssignOri(xloc,yloc,zloc,RotAngle):
    import section
    import regionToolset
    import displayGroupMdbToolset as dgm
    import part
    import material
    import assembly
    import step
    import interaction
    import load
    import mesh
    import optimization
    import job
    import sketch
    import visualization
    import xyPlot
    import displayGroupOdbToolset as dgo
    import connectorBehavior
    p = mdb.models['Model-1'].parts['Part-1']
    c = p.cells
    cells = c.findAt(((xloc, yloc, zloc), ))
    region = regionToolset.Region(cells=cells)
    orientation=None
    p.MaterialOrientation(region=region, orientationType=SYSTEM, 
		axis=AXIS_3, localCsys=orientation, fieldName='', 
		additionalRotationType=ROTATION_ANGLE, additionalRotationField='',
		angle=RotAngle, stackDirection=STACK_3)


currentDir = os.path.abspath(os.getcwd())
caeName = 'LVI_1A'
stepType = 'Explicit'
mdb.saveAs(pathName=currentDir + '/' + caeName + '.cae')
#
#
#
# put laminate configuration below
NoPly = 16
NoIntf = 15
PlyThickness = 0.2375
IntfThickness = 0.00001
TotWidth = 100.0
TotLength = 150.0
TotTime = 200.0
TotThickness = float(NoPly)*PlyThickness+float(NoIntf)*IntfThickness
OrientationList = [45.0,90.0,-45.0,0.0,45.0,90.0,-45.0,0.0,0.0,-45.0,90.0,45.0,0.0,-45.0,90.0,45.0]
impactVelocity = -2049

#python scripting starts here
session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)
myModel1 = mdb.models['Model-1'] 
myModel1.ConstrainedSketch(name='__profile__', sheetSize=200.0)
myModel1.sketches['__profile__'].rectangle(point1=(0.0, 0.0), 
    point2=(TotLength, TotWidth))
myModel1.Part(dimensionality=THREE_D, name='Part-1', type=
    DEFORMABLE_BODY)
myModel1.parts['Part-1'].BaseSolidExtrude(depth=TotThickness, sketch=
    myModel1.sketches['__profile__'])
del myModel1.sketches['__profile__']
# Save by pulungd on 2018_10_15-10.06.33; build 6.14-5 2015_08_18-17.37.49 135153
PlyCounter = 0
IntfCounter = 0
for Ori in range(len(OrientationList)-1):
	PlyCounter = PlyCounter + 1
	MyOffset = float(PlyCounter)*PlyThickness+float(IntfCounter)*IntfThickness
	MakeDatumbyOffset(MyOffset,TotLength/2,TotWidth/2,0.0)
	print 'Ply', OrientationList[Ori], 'at offset', MyOffset
	if (OrientationList[PlyCounter] != OrientationList[PlyCounter-1]):
		IntfCounter = IntfCounter + 1
		MakeDatumbyOffset(MyOffset+IntfThickness,TotLength/2,TotWidth/2,0.0)
		print 'Interface', IntfCounter, 'at offset', MyOffset+IntfThickness
		


mdb.save()
p1 = myModel1.parts['Part-1']
session.viewports['Viewport: 1'].setValues(displayedObject=p1)
session.viewports['Viewport: 1'].view.setValues(session.views['Bottom'])
DatumNo = 1
for PartitionNo in range(PlyCounter+IntfCounter):
	DatumNo = DatumNo + 1
	zheight = TotThickness-PlyThickness/2.0
	PartitionIt(zheight,DatumNo)
	print PartitionNo+1,'-st Datum Partition'


# Making tools, pin holder, base plate fixture, impactor
BallRad = 8.0
PinRad = 7.0
pinHeight = 4.0

createBaseFixture(100.0,150.0,75.0,125.0,'Part-2-BasePlate')
createPinHolder(PinRad,pinHeight,'Part-3-PinHolder')
createImpactorBall(BallRad,'Part-4-Impactor')

p2 = myModel1.parts['Part-2-BasePlate']
p3 = myModel1.parts['Part-3-PinHolder']
p4 = myModel1.parts['Part-4-Impactor']

# Making materials, section and assign it
# making ply material
mat1 = myModel1.Material(name='Ply_matl')
mat1.Density(table=((1.59e-9, ), ))
mat1.Elastic(type=ENGINEERING_CONSTANTS, table=((131610.0, 9238.0, 9238.0, 0.302, 
    0.302, 0.35, 4826.0, 4826.0, 3548.0), ))
mat1.HashinDamageInitiation(table=(( 2063.0, 1484.0, 63.0, 267.0, 91.0, 133.0), ))
mat1.hashinDamageInitiation.DamageEvolution( table=((106.3, 81.5, 0.28, 0.79), ), type=ENERGY)
mat1.hashinDamageInitiation.DamageStabilization(
    fiberCompressiveCoeff=0.0001, fiberTensileCoeff=0.0001, 
    matrixCompressiveCoeff=0.0001, matrixTensileCoeff=0.0001)
myModel1.HomogeneousShellSection(idealization=NO_IDEALIZATION, 
    integrationRule=SIMPSON, material='Ply_matl', name='Ply_section', 
    nodalThicknessField='', numIntPts=5, poissonDefinition=DEFAULT, 
    preIntegrate=OFF, temperature=GRADIENT, thickness=PlyThickness, thicknessField=''
    , thicknessModulus=None, thicknessType=UNIFORM, useDensity=OFF)

#making cohesive interface material
mat2 = myModel1.Material(name='Intf_matl')
mat2.Density(table=((0.1111, ), ))
mat2.Elastic(type=TRACTION, table=((
    50000.0, 50000.0, 50000.0), ))
mat2.QuadsDamageInitiation(table=((
    26.26, 31.89, 31.89), ))
mat2.quadsDamageInitiation.DamageEvolution(
    type=ENERGY, mixedModeBehavior=BK, power=1.45, table=((0.28, 0.79, 
    0.79), ))
mat2.quadsDamageInitiation.DamageStabilizationCohesive(
    cohesiveCoeff=0.0001*TotTime)
myModel1.CohesiveSection(name='Intf_section', 
    material='Intf_matl', response=TRACTION_SEPARATION, 
    outOfPlaneThickness=None)
mdb.save()
#
# Assigning ply material section
cellsPly = []
cellsIntf = []
c1 = p1.cells
for PlyNo in range(len(OrientationList)):
	xloc = TotLength/2.0
	yloc = TotWidth/2.0
	zloc = float(PlyNo)*PlyThickness+PlyThickness/2.0
	mySetName = 'Ply_'+str(PlyNo+1)+'_set'
	mySection = 'Ply_section'
	AssignMatl(xloc,yloc,zloc,mySetName,mySection)
	RotAngle = OrientationList[PlyNo]
	AssignOri(xloc,yloc,zloc,RotAngle)
	print mySetName,' has been assigned section', mySection
	cellsPly.append(c1.findAt(((xloc, yloc, zloc), )))


#
# Assigning interface material section
PlyCounter = 0
IntfCounter = 0
xloc = TotLength/2.0
yloc = TotWidth/2.0
for Ori in range(len(OrientationList)-1):
	PlyCounter = PlyCounter + 1
	MyOffset = float(PlyCounter)*PlyThickness+float(IntfCounter)*IntfThickness
	if (OrientationList[PlyCounter] != OrientationList[PlyCounter-1]):
		IntfCounter = IntfCounter + 1
		zloc = MyOffset+IntfThickness/2.0
		mySetName = 'Intf_'+str(IntfCounter)+'_set'
		mySection = 'Intf_section'
		AssignMatl(xloc,yloc,zloc,mySetName,mySection)
        cellsIntf.append(c1.findAt(((xloc, yloc, zloc), )))
        print zloc
        print mySetName,' has been assigned section', mySection


myModel1.HomogeneousShellSection(idealization=NO_IDEALIZATION, 
    integrationRule=SIMPSON, material='Ply_matl', name='Section-3', 
    nodalThicknessField='', numIntPts=5, poissonDefinition=DEFAULT, 
    preIntegrate=OFF, temperature=GRADIENT, thickness=0.2375, thicknessField=''
    , thicknessModulus=None, thicknessType=UNIFORM, useDensity=OFF)

# Making assembly/instance of the part
myModel1.rootAssembly.DatumCsysByDefault(CARTESIAN)
myAssembly = myModel1.rootAssembly

#creating the laminate instance
myAssembly.Instance(dependent=ON, name='Part-1-1', 
    part=p1)

#creating the impactor instance
myAssembly.Instance(dependent=ON, name=
    'Part-4-Impactor-1', part=p4)
myAssembly.translate(instanceList=('Part-4-Impactor-1', 
    ), vector=(TotLength/2.0, TotWidth/2.0, NoPly*PlyThickness+1.0+BallRad))

#creating the base plate
myAssembly.Instance(dependent=ON, name=
    'Part-2-BasePlate-1', part=p2)
myAssembly.rotate(angle=-90.0, axisDirection=(0.0, 0.0, 
    1), axisPoint=(0.0, 0.0, 0.0), instanceList=('Part-2-BasePlate-1', ))
myAssembly.translate(instanceList=(
    'Part-2-BasePlate-1', ), vector=(0.0, TotWidth, 0.0))

#creating 4 pin holders
myAssembly.Instance(dependent=ON, name=
    'Part-3-PinHolder-1', part=p3)
myAssembly.instances['Part-3-PinHolder-1'].translate(
    vector=(1.0/6.0*TotLength, PinRad, NoPly*PlyThickness))
myAssembly.Instance(dependent=ON, name=
    'Part-3-PinHolder-2', part=p3)
myAssembly.instances['Part-3-PinHolder-2'].translate(
    vector=(5.0/6.0*TotLength, PinRad, NoPly*PlyThickness))
myAssembly.Instance(dependent=ON, name=
    'Part-3-PinHolder-3', part=p3)
myAssembly.instances['Part-3-PinHolder-3'].translate(
    vector=(1.0/6.0*TotLength, TotWidth - PinRad, NoPly*PlyThickness))
myAssembly.Instance(dependent=ON, name=
    'Part-3-PinHolder-4', part=p3)
myAssembly.instances['Part-3-PinHolder-4'].translate(
    vector=(5.0/6.0*TotLength,TotWidth - PinRad, NoPly*PlyThickness))



#setting up step, field and history output #dynamic #explicit case
myModel1.ExplicitDynamicsStep(name='Step-1', previous='Initial', 
    timePeriod=0.005)
myModel1.fieldOutputRequests['F-Output-1'].setValues(numIntervals=
    200, rebar=INCLUDE, sectionPoints=DEFAULT, variables=('S', 'SVAVG', 'PE', 
    'PEVAVG', 'PEEQ', 'PEEQVAVG', 'LE', 'U', 'V', 'A', 'RF', 'CSTRESS', 
    'DAMAGEC', 'DAMAGET', 'DAMAGEFT', 'DAMAGEFC', 'DAMAGEMT', 'DAMAGEMC', 
    'DAMAGESHR', 'SDEG', 'DMICRT', 'EVF', 'IVOL', 'STATUS'))
myModel1.historyOutputRequests['H-Output-1'].setValues(variables=(
    'ALLAE', 'ALLCD', 'ALLDMD', 'ALLFD', 'ALLIE', 'ALLKE', 'ALLPD', 'ALLSE', 
    'ALLVD', 'ALLWK', 'ALLMW', 'ETOTAL'))

# myAssembly.Surface(name='m_Surf-2', side1Edges=
#     myAssembly.instances['ply0_bot_1'].edges.getByBoundingBox(
#     	-0.1,plythickness/2.0-0.1,0,length_x+0.1,plythickness/2.0+0.1,0.0))

# myAssembly.Set(faces=
#     myAssembly.instances['Part-1-1'].faces.findAt(((
#     0.0, 6.6, 3.825097), ), ((0.0, 13.199999, 3.9101), ), ((0.0, 6.6, 
#     3.570087), ), ((0.0, 6.6, 3.74009), ), ((0.0, 6.6, 3.23008), ), ((0.0, 6.6, 
#     3.48508), ), ((0.0, 6.6, 2.97507), ), ((0.0, 6.6, 3.060077), ), ((0.0, 6.6, 
#     2.72006), ), ((0.0, 6.6, 2.805067), ), ((0.0, 6.6, 2.46505), ), ((0.0, 6.6, 
#     2.550057), ), ((0.0, 6.6, 1.95505), ), ((0.0, 6.6, 2.21005), ), ((0.0, 6.6, 
#     1.530047), ), ((0.0, 6.6, 1.70005), ), ((0.0, 6.6, 1.275037), ), ((0.0, 
#     6.6, 1.44504), ), ((0.0, 6.6, 1.020027), ), ((0.0, 6.6, 1.19003), ), ((0.0, 
#     6.6, 0.68002), ), ((0.0, 6.6, 0.93502), ), ((0.0, 6.6, 0.42501), ), ((0.0, 
#     6.6, 0.510017), ), ((0.0, 6.6, 0.17), ), ((0.0, 6.6, 0.255007), ), ), name=
#     'Set-1')


# creating encastre boundary conditions for pins reference points
myAssembly.Set(name='Set-Fix-RP_Pins', referencePoints=
    ((myAssembly.instances['Part-3-PinHolder-3'].referencePoints[3], 
    ), (myAssembly.instances['Part-3-PinHolder-4'].referencePoints[3], 
    ), (myAssembly.instances['Part-3-PinHolder-2'].referencePoints[3], 
    ), (myAssembly.instances['Part-3-PinHolder-1'].referencePoints[3], 
    )))
myModel1.EncastreBC(createStepName='Initial', localCsys=None, name=
    'BC-1', region=myAssembly.sets['Set-Fix-RP_Pins'])

# creating pin boundary conditions for the impactor, except for z direction.
myAssembly.Set(name='Set-RP-Impactor', referencePoints=
    (myAssembly.instances['Part-4-Impactor-1'].referencePoints[2], 
    ))
myModel1.DisplacementBC(amplitude=UNSET, createStepName='Initial', 
    distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
    'BC-2', region=mdb.models['Model-1'].rootAssembly.sets['Set-RP-Impactor'], 
    u1=0.0, u2=0.0, u3=UNSET, ur1=0.0, ur2=0.0, ur3=0.0)

mdb.models['Model-1'].Velocity(distributionType=MAGNITUDE, field='', name=
    'Predefined Field-1', omega=0.0, region=
    mdb.models['Model-1'].rootAssembly.sets['Set-RP-Impactor'], velocity1=0.0, 
    velocity2=0.0, velocity3=-impactVelocity)


# creating BC for the baseplate
myAssembly.Set(name='Set-RP-Baseplate', 
    referencePoints=(myAssembly.instances['Part-2-BasePlate-1'].referencePoints[2], ))
myModel1.EncastreBC(createStepName='Initial', localCsys=None, 
    name='BC-3', region=myAssembly.sets['Set-RP-Baseplate'])


	
# generate meshing for each part 	Part 1
p1.seedPart(deviationFactor=0.1, minSizeFactor=0.1, size=2.5)
p1.generateMesh()
myAssembly.regenerate()

#set element types to continuum shell elements
p1.setElementType(elemTypes=(ElemType(
    elemCode=SC8R, elemLibrary=EXPLICIT, secondOrderAccuracy=ON, 
    hourglassControl=ENHANCED, elemDeletion=ON, maxDegradation=0.999), 
    ElemType(elemCode=SC6R, elemLibrary=EXPLICIT, secondOrderAccuracy=ON), 
    ElemType(elemCode=UNKNOWN_TET, elemLibrary=EXPLICIT)), regions=(cellsPly))


#set element types for cohesive elements
p1.setElementType(elemTypes=(ElemType(
    elemCode=COH3D8, elemLibrary=EXPLICIT, elemDeletion=ON, 
    maxDegradation=0.999), ElemType(elemCode=COH3D6, elemLibrary=EXPLICIT, 
    elemDeletion=ON, maxDegradation=0.999), ElemType(elemCode=UNKNOWN_TET, 
    elemLibrary=EXPLICIT)), regions=(cellsIntf))


#
# Assigning stacking direction for ply elements
PlyCounter = 0
IntfCounter = 0
xloc = TotLength/2.0
yloc = TotWidth/2.0
zref = TotThickness
for Ori in range(len(OrientationList)): #note the counter here diferent with interface one
    MyOffset = float(PlyCounter)*PlyThickness+float(IntfCounter)*IntfThickness
    PlyCounter = PlyCounter + 1
    zloc = MyOffset+PlyThickness/2.0
    AssignStackDir(xloc,yloc,zloc,xloc,yloc,zref)
    if PlyCounter != len(OrientationList): #this is to avoid out of index for the last ply
        if (OrientationList[PlyCounter] != OrientationList[PlyCounter-1]): #differs also with interface one
            IntfCounter = IntfCounter + 1
    print 'Success assigning stack direction for ply ', PlyCounter
    


#
# Assigning stacking direction for interface elements
PlyCounter = 0
IntfCounter = 0
xloc = TotLength/2.0
yloc = TotWidth/2.0
zref = TotThickness
for Ori in range(len(OrientationList)-1):
	PlyCounter = PlyCounter + 1
	MyOffset = float(PlyCounter)*PlyThickness+float(IntfCounter)*IntfThickness
	if (OrientationList[PlyCounter] != OrientationList[PlyCounter-1]):
		IntfCounter = IntfCounter + 1
		zloc = MyOffset+IntfThickness/2.0
		AssignStackDir(xloc,yloc,zloc,xloc,yloc,zref)
		print 'Success assigning stack direction for interface ', IntfCounter


# meshing part 2
p2.seedPart(deviationFactor=0.1, minSizeFactor=0.1, size=2.5)
p2.setElementType(elemTypes=(ElemType(elemCode=R3D4, elemLibrary=EXPLICIT), 
    ElemType(elemCode=R3D3, elemLibrary=EXPLICIT)), regions=(
    p2.faces.findAt(((91.666667, 50.0, 0.0), )), ))
p2.setMeshControls(elemShape=QUAD, regions=p2.faces.findAt(((
    91.666667, 50.0, 0.0), )))
p2.generateMesh()

# meshing part 3
p3.seedPart(deviationFactor=0.1, minSizeFactor=0.1, size=1.0)
p3.setElementType(elemTypes=(ElemType(elemCode=R3D4, elemLibrary=EXPLICIT), 
    ElemType(elemCode=R3D3, elemLibrary=EXPLICIT)), regions=(
    p3.faces.findAt(((-6.982941, -0.488394, 1.333333), ), 
    ((0.0, -6.898022, 4.0), ), ((0.0, -6.898022, 0.0), ), ), ))
p3.generateMesh()


# meshing part 4
p4.seedPart(deviationFactor=0.1, minSizeFactor=0.1, size=0.5)
p4.setElementType(elemTypes=(ElemType(elemCode=R3D4, elemLibrary=EXPLICIT), 
    ElemType(elemCode=R3D3, elemLibrary=EXPLICIT)), regions=(
    p4.faces.findAt(((3.202587, -7.301164, -0.660635), )), ))
p4.generateMesh()



# to generate mass of the reference point
p4.engineeringFeatures.PointMassInertia( alpha=0.0, composite=0.0, 
    mass=0.005, name='Inertia-1', region= p4.sets['Set-RP-Part-4-Impactor'])


myAssembly.regenerate()
myModel1.ContactProperty('IntProp-1')
myModel1.interactionProperties['IntProp-1'].TangentialBehavior(
    dependencies=0, directionality=ISOTROPIC, elasticSlipStiffness=None, 
    formulation=PENALTY, fraction=0.005, maximumElasticSlip=FRACTION, 
    pressureDependency=OFF, shearStressLimit=None, slipRateDependency=OFF, 
    table=((0.3, ), ), temperatureDependency=OFF)
myModel1.interactionProperties['IntProp-1'].NormalBehavior(
    allowSeparation=ON, constraintEnforcementMethod=DEFAULT, 
    pressureOverclosure=HARD)
myAssembly.Surface(name='m_Surf-1', side1Faces=
    myAssembly.instances['Part-4-Impactor-1'].faces.findAt(
    ((78.202587, 42.698836, 12.139365), )))
myAssembly.Set(faces=
    myAssembly.instances['Part-1-1'].faces.findAt(((
    50.0, 33.333333, 3.80015), )), name='s_Set-4')
myModel1.SurfaceToSurfaceContactExp(clearanceRegion=None, 
    createStepName='Initial', datumAxis=None, initialClearance=OMIT, 
    interactionProperty='IntProp-1', master=
    myAssembly.surfaces['m_Surf-1'], 
    mechanicalConstraint=PENALTY, name='Int-1', slave=
    myAssembly.sets['s_Set-4'], sliding=FINITE)
myAssembly.Surface(name='m_Surf-2', side1Faces=
    myAssembly.instances['Part-2-BasePlate-1'].faces.findAt(
    ((50.0, 8.333333, 0.0), )))
myAssembly.Surface(name='s_Surf-2', side1Faces=
    myAssembly.instances['Part-1-1'].faces.findAt(((
    100.0, 33.333333, 0.0), )))
myModel1.SurfaceToSurfaceContactExp(clearanceRegion=None, 
    createStepName='Initial', datumAxis=None, initialClearance=OMIT, 
    interactionProperty='IntProp-1', master=
    myAssembly.surfaces['m_Surf-2'], 
    mechanicalConstraint=KINEMATIC, name='Int-2', slave=
    myAssembly.surfaces['s_Surf-2'], sliding=FINITE)
# Save by pulungd on 2021_02_10-22.58.23; build 2017 2016_09_28-04.54.59 126836


# #left support pin
# MakeDatumbyOffset(12.0,0.0,TotWidth/2.0,TotThickness/2.0)
# MakeDatumbyOffset(17.0,0.0,TotWidth/2.0,TotThickness/2.0)
# MakeDatumbyOffset(22.0,0.0,TotWidth/2.0,TotThickness/2.0)
# #left loading pin
# MakeDatumbyOffset(31.0,0.0,TotWidth/2.0,TotThickness/2.0)
# MakeDatumbyOffset(36.0,0.0,TotWidth/2.0,TotThickness/2.0)
# MakeDatumbyOffset(41.0,0.0,TotWidth/2.0,TotThickness/2.0)
# #x-symmetry datum plane
# MakeDatumbyOffset(55.0,0.0,TotWidth/2.0,TotThickness/2.0)
# #right loading pin
# MakeDatumbyOffset(69.0,0.0,TotWidth/2.0,TotThickness/2.0)
# MakeDatumbyOffset(74.0,0.0,TotWidth/2.0,TotThickness/2.0)
# MakeDatumbyOffset(79.0,0.0,TotWidth/2.0,TotThickness/2.0)
# #right support pin
# MakeDatumbyOffset(88.0,0.0,TotWidth/2.0,TotThickness/2.0)
# MakeDatumbyOffset(93.0,0.0,TotWidth/2.0,TotThickness/2.0)
# MakeDatumbyOffset(98.0,0.0,TotWidth/2.0,TotThickness/2.0)

	

mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
    explicitPrecision=DOUBLE_PLUS_PACK, getMemoryFromAnalysis=True, historyPrint=OFF, 
    memory=90, memoryUnits=PERCENTAGE, model='Model-1', modelPrint=OFF, 
    multiprocessingMode=DEFAULT, name='LVI_1A', nodalOutputPrecision=FULL, 
    numCpus=6, numDomains=6, numGPUs=0, queue=None, resultsFormat=ODB, scratch=
    '', type=ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)

#mdb.jobs['Job-1'].submit(consistencyChecking=OFF)