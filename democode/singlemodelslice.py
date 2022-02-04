import sys
sys.path.append(r"/home/goatchurch/geom3d/barmesh")
from basicgeo import P2, P3, Partition1, Along
import barmesh, implicitareaballoffset, implicitareacyloffset
from mainfunctions import BarMeshContoursF, PlotBarmesh, NestContours
from barmeshslicer import BarMeshSlicer
from tribarmes import TriangleBarMesh

sendactivity("clearalltriangles")
sendactivity("clearallcontours")
sendactivity("clearallpoints")

fname = "/home/goatchurch/geom3d/barmesh/stlsamples/frameguide.stl"
tbm = TriangleBarMesh(fname)
iaoffset = implicitareaballoffset.ImplicitAreaBallOffset(tbm)
sendactivity("triangles", codetriangles=tbm.GetBarMeshTriangles())

rad = 2.5
rex = rad + 2.5
xpart = Partition1(tbm.xlo-rex, tbm.xhi+rex, 45)
ypart = Partition1(tbm.ylo-rex, tbm.yhi+rex, 37)
zlevel = Along(0.5, tbm.zlo, tbm.zhi)

bm = barmesh.BarMesh()
bm.BuildRectBarMesh(xpart, ypart, zlevel)
rd2 = max(xpart.vs[1]-xpart.vs[0], ypart.vs[1]-ypart.vs[0], rad*1.5) + 0.1
bms = BarMeshSlicer(bm, iaoffset, rd=rad, rd2=rd2, contourdotdiff=0.95, contourdelta=0.05, lamendgap=0.001)
bms.fullmakeslice()

conts, topbars = BarMeshContoursF(bm, barmesh.PZ_BEYOND_R)
sendactivity(contours=conts)

PlotBarmesh(bm, sendactivity)


# now demonstrate that we can nest the contours and plot the inner and outer ones
contnest = NestContours(topbars, barmesh.PZ_BEYOND_R)
mconts = dict((topbar.midcontournumber, cont)  for cont, topbar in zip(conts, topbars))

sendactivity("contours", contours=[mconts[cn]  for cn, (izone, outxn, innlist) in contnest.items()  if izone == barmesh.PZ_BEYOND_R], materialnumber=0)
sendactivity("contours", contours=[mconts[cn]  for cn, (izone, outxn, innlist) in contnest.items()  if izone != barmesh.PZ_BEYOND_R], materialnumber=2)

