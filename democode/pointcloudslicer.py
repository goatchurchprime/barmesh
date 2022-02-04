from basicgeo import P2, P3, Partition1, Along
import barmesh, implicitareaballoffset, implicitareacyloffset
from mainfunctions import BarMeshContoursF, PlotBarmesh, NestContours
from barmeshslicer import BarMeshSlicer
from tribarmes import TriangleBarMesh
from tribarmes import MakeTriangleBoxing
import struct

def GetTBMfromPoints():
    fname = "/home/goatchurch/caving/Freeze.ply"
    fin = open(fname, "rb")
    k = fin.read()
    eh = b"end_header\n"
    meh = k.find(eh) + len(eh)
    k1 = k[meh:]
    n = 1320559
    tbm = TriangleBarMesh()
    for i in range(n):
        tbm.NewNode(P3(*struct.unpack("<fff", k1[i*12:i*12+12])))
    return tbm
tbm = GetTBMfromPoints()
print("Made Triangle bar mesh")
iaoffset = implicitareaballoffset.ImplicitAreaBallOffset(tbm)
print("Made implicitarea offset")
#sendactivity(points=[struct.unpack("<fff", k1[i*12:i*12+12])  for i in range(n)])

def GetSlice(zlevel):
    rad = 1.1
    rex = rad + 2.5
    xpart = Partition1(tbm.xlo-rex, tbm.xhi+rex, 45)
    ypart = Partition1(tbm.ylo-rex, tbm.yhi+rex, 37)

    bm = barmesh.BarMesh()
    bm.BuildRectBarMesh(xpart, ypart, zlevel)
    rd2 = max(xpart.vs[1]-xpart.vs[0], ypart.vs[1]-ypart.vs[0], rad*1.5) + 0.1
    bms = BarMeshSlicer(bm, iaoffset, rd=rad, rd2=rd2, contourdotdiff=0.95, contourdelta=0.05, lamendgap=0.001)
    bms.fullmakeslice()

    conts, topbars = BarMeshContoursF(bm, barmesh.PZ_BEYOND_R)
    print("z:%.3f pztime:%.3f/%d cbtime:%.3f(+%.3f)/%d remains: %3f" % (zlevel, bms.pztime, bms.pzcalls, bms.cbtime, bms.cbztime, bms.cbcalls, (bms.totalproctime - bms.pztime - bms.cbtime)))
    return conts, topbars

zlevels = list(range(int(tbm.zlo+1), int(tbm.zhi+1)))
outercontslist = [ ]
for zlevel in zlevels[0:]:
    print("doing level", zlevel)
    conts, topbars = GetSlice(zlevel)
    contnest = NestContours(topbars, barmesh.PZ_BEYOND_R)
    mconts = dict((topbar.midcontournumber, cont)  for cont, topbar in zip(conts, topbars))
    outerconts = [mconts[cn]  for cn, (izone, outxn, innlist) in contnest.items()  if izone == barmesh.PZ_BEYOND_R and outxn == -1]
    innerconts = [mconts[cn]  for cn, (izone, outxn, innlist) in contnest.items()  if not (izone == barmesh.PZ_BEYOND_R and outxn == -1)]
    outercontslist.append((zlevel, outerconts))
#    sendactivity(contours=outerconts)
#    sendactivity(contours=innerconts, materialnumber=2)

import pickle
pickle.dump(outercontslist, open("/home/goatchurch/caving/Freeze.conts", "wb"))
    
# export PYTHONPATH=/home/goatchurch/geom3d/barmesh
# pypy pointcloudslicer.py
