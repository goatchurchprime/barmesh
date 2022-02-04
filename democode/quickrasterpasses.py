import sys
from basicgeo import P2, P3, Partition1, Along
import barmesh, implicitareaballoffset, implicitareacyloffset
import mainfunctions
from barmeshslicer import BarMeshSlicer

from tribarmes import TriangleBarMesh, MakeTriangleBoxing

fname = "/home/goatchurch/geom3d/barmesh/stlsamples/frameguide.stl"
tbm = TriangleBarMesh(fname)
ibo = implicitareaballoffset.ImplicitAreaBallOffset(tbm)
sendactivity("triangles", codetriangles=tbm.GetBarMeshTriangles())

rad = 2.5
rex = rad*2
xpart = Partition1(tbm.xlo-rex, tbm.xhi+rex, 27)
ypart = Partition1(tbm.ylo-rex, tbm.yhi+rex, 96)

for x in xpart.vs:
    cont = [ ]
    for y in ypart.vs:
        p0 = P3(x, y, tbm.zhi + rad*2)
        vp = P3(12.1, 12.3, tbm.zlo - tbm.zhi - rad*4)
        lam = ibo.Cutpos(p0, vp, None, rad)
        cont.append(p0 + vp*min(1,lam))
    sendactivity(contours=[cont])
