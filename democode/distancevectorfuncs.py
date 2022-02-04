import sys, imp
sys.path.append(r"/home/goatchurch/geom3d/barmesh")
from basicgeo import P2, P3, Partition1, Along
import barmesh, triangleboxing, mainfunctions
import implicitareaballoffset, implicitareacyloffset

sendactivity("clearalltriangles")
fname = "/home/goatchurch/geom3d/stlfiles/impellor1.stl"
#tbm = mainfunctions.TBMfromSTL(fname, (lambda t: (t[0], -t[2], t[1])))
#z = 6.5
fname = "/home/goatchurch/geom3d/stlfiles/barobox.stl"
tbm = mainfunctions.TBMfromSTL(fname)

sendactivity("triangles", codetriangles=tbm.GetBarMeshTriangles())
ibo = implicitareaballoffset.ImplicitAreaBallOffset(tbm)

ps = [ ]
conts = [ ]
for i in range(20):
    for j in range(20):
        p = P3(Along(i/19, tbm.xlo, tbm.xhi), Along(j/19, tbm.ylo, tbm.yhi), tbm.zhi+5)
        pz = barmesh.PointZone(0, 30, None)
        ibo.DistP(pz, p)
        ps.append(p)
        if pz.v:
            conts.append([p, p+pz.v*0.9])
        
sendactivity(points=ps)
sendactivity(contours=conts)

# find the vector values too
conts = [ ]
for p in ps:
    vp = P3(6,1,-20)
    lam = ibo.Cutpos(p, vp, None, 3)
    if lam != 2:
        conts.append([p, p+vp*lam])
sendactivity(contours=conts, materialnumber=1)


