from basicgeo import P2, P3, Partition1, Along
import barmesh, implicitareaballoffset, implicitareacyloffset
from mainfunctions import BarMeshContoursF, PlotBarmesh, NestContours
from barmeshslicer import BarMeshSlicer
from tribarmes import TriangleBarMesh
from tribarmes.trianglebarmesh import TriangleBar
from tribarmes import MakeTriangleBoxing
import struct

import pickle
outercontslist = pickle.load(open("/home/goatchurch/caving/Freeze.conts", "rb"))
for z, conts in outercontslist:
    sendactivity(contours=conts)


conts = [[(0,0,0), (1,2,0), (3,1,0), (0,0,0)]]
def GetTBMfromConts(conts):
    tbm = TriangleBarMesh()
    for cont in conts:
        i0 = len(tbm.nodes)
        for p in cont[:-1]:
            tbm.NewNode(P3(p[0], p[1], p[2]))
        i1 = len(tbm.nodes)
        for i in range(i0+1, i1):
            tbm.bars.append(TriangleBar(tbm.nodes[i-1], tbm.nodes[i]))
        tbm.bars.append(TriangleBar(tbm.nodes[i0], tbm.nodes[i1-1]))
    return tbm
tbm = GetTBMfromConts(conts)
tbm.zlo, tbm.zhi
iaoffset = implicitareaballoffset.ImplicitAreaBallOffset(tbm)

def GetSlice(tbm):
    rad = 1.1+0.02
    rex = rad + 2.5
    xpart = Partition1(tbm.xlo-rex, tbm.xhi+rex, 45)
    ypart = Partition1(tbm.ylo-rex, tbm.yhi+rex, 37)

    bm = barmesh.BarMesh()
    bm.BuildRectBarMesh(xpart, ypart, tbm.zlo)
    rd2 = max(xpart.vs[1]-xpart.vs[0], ypart.vs[1]-ypart.vs[0], rad*1.5) + 0.1
    bms = BarMeshSlicer(bm, iaoffset, rd=rad, rd2=rd2, contourdotdiff=0.95, contourdelta=0.05, lamendgap=0.001)
    bms.fullmakeslice()

    conts, topbars = BarMeshContoursF(bm, barmesh.PZ_BEYOND_R)
    print("z:%.3f pztime:%.3f/%d cbtime:%.3f(+%.3f)/%d remains: %3f" % (zlevel, bms.pztime, bms.pzcalls, bms.cbtime, bms.cbztime, bms.cbcalls, (bms.totalproctime - bms.pztime - bms.cbtime)))
    return conts, topbars
