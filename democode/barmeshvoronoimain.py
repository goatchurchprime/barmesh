from basicgeo import P2, Partition1, I1, Along
from vor2d import VBconstructionbarmesh
from vor2d import VBdrivegeometry, VBconstructiongeometry
from vor2d import VBdiagramcreation, InsertColPointNodes
from vor2d import PlotVBdiagram

#import voronoiexpt, voronoisegment
#voronoiexpt.sendactivity = sendactivity
#voronoisegment.sendactivity = sendactivity

sendactivity("clearalltriangles")
sendactivity("clearallcontours")
sendactivity("clearallpoints")

pts = [P2(0,0), P2(10,-2), P2(9,5), P2(5.4,4.3), P2(2,6)]
pts1 = [P2(3.9,3.12),P2(7.33,1.28),P2(4.49,0.63)]
#pts1.reverse()

e = 0.01
#pts = [ P2(*p)  for p in [(-0.5,4.7),(1.8,7.0),(3.2,7.0+e),(4.7,7.0),(5.3,6.1),(7.6,3.8),(8.6,3.8),(10.0,3.6),(10.7,2.6),(10.6,0.4),(9.4,-0.5),(7.2,-0.6),(4.9,0.1),(3.7,1.0),(3.6,0.0),(2.6,-1.4),(1.1,-1.3),(-0.2,-0.5),(-0.6,0.8),(0.6,2.2),(2.0,2.7)] ]
pts.append(pts[0])
pts1.append(pts1[0])

conts = [pts, pts1]
urg = I1.AbsorbList(p.u  for cont in conts  for p in cont)
vrg = I1.AbsorbList(p.v  for cont in conts  for p in cont)

upart = Partition1(urg.lo-5.01, urg.hi+5, 3+1)
vpart = Partition1(vrg.lo-5.03, vrg.hi+5, 3+3)

#xpart = Partition1(-5, 15, 30+1)
#ypart = Partition1(-5, 10, 30+3)

bm = VBconstructionbarmesh()
bm.BuildRectBarMesh(xpart, ypart)

vbs = VBdrivegeometry(I1(urg.lo-50, urg.hi+50), I1(vrg.lo-50, vrg.hi+50))
[node.p  for node in bm.nodes]



sendactivity(contours=[pts], materialnumber=2)
sendactivity(points=pts)
vbs.AddPoly(pts)

sendactivity(contours=[pts1], materialnumber=2)
vbs.AddPoly(pts1)
vbs.AddPolyDone(I1(upart.lo-1, upart.hi+1), I1(vpart.lo-1, vpart.hi+1))

vbs.sendactivity = sendactivity

vbc = VBconstructiongeometry(vbs, bm)
vbc.goconstruct(maxDnextvbarstack=1600)
print("Construct time %f vbscloset %f%% calls=%d" % (vbc.vbstotaltime, vbc.vbsclosesttime*100/vbc.vbstotaltime, vbc.vbsclosestcalls))

#vbc.PlotBarMesh(sendactivity, False)

vbd = VBdiagramcreation(vbs, vbc)
vbd.BuildVBDC()
vbm = vbd.vbm

sendactivity(contours=PlotVBdiagram(vbs, vbm, rc=-1, coslimit=1))

Nn = len(vbm.nodes)
InsertColPointNodes(vbs, vbm)
sendactivity(points=[node.p  for node in vbm.nodes[Nn:]], materialnumber=3)

