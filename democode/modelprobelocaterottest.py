import sys, random, math, time
sys.path.append(r"/home/goatchurch/geom3d/barmesh")
from basicgeo import P2, P3, Partition1, Along, Quat
import barmesh, triangleboxing, mainfunctions
import implicitareaballoffset, implicitareacyloffset
import scipy.optimize

sendactivity("clearalltriangles")
sendactivity("clearallpoints")
sendactivity("clearallcontours")

fname = "/home/goatchurch/geom3d/barmesh/stlsamples/frameguide.stl"
tbm = mainfunctions.TBMfromSTL(fname)
sendactivity("triangles", codetriangles=tbm.GetBarMeshTriangles())

ibo = implicitareaballoffset.ImplicitAreaBallOffset(tbm)

proberad = 2
Nprobes = 10

    
def PointsInRadRange(midpoint, r0, r1, n):
    res = [ ]
    while len(res) < n:
        p = P3(random.uniform(midpoint.x-r1, midpoint.x+r1), random.uniform(midpoint.y-r1, midpoint.y+r1), random.uniform(midpoint.z-r1, midpoint.z+r1))
        if r0 < (p - midpoint).Len() < r1:  res.append(p)
    return res
    
def MakeProbePoints():
    midpoint = P3((tbm.xlo+tbm.xhi)/2, (tbm.ylo+tbm.yhi)/2, (tbm.zlo+tbm.zhi)/2)
    modelrad = (P3(tbm.xhi, tbm.yhi, tbm.zhi) - midpoint).Len()

    probefrompoints = PointsInRadRange(midpoint, modelrad+5, modelrad*3, Nprobes)
    probetopoints = PointsInRadRange(midpoint, 0, modelrad/2, Nprobes)
    conts = [ ]
    for probefrompoint, probetopoint in zip(probefrompoints, probetopoints):
        vp = probetopoint - probefrompoint
        lam = ibo.Cutpos(probefrompoint, vp, None, proberad)
        if lam != 2:
            conts.append([probefrompoint, probefrompoint+vp*lam])
    sendactivity(contours=conts, materialnumber=1)
    probepoints = [p  for fp, p in conts]
    sendactivity(points=probepoints)
    return probepoints

def MeasureDistanceToModel(p):
    pz = barmesh.PointZone(0, 30, None)
    ibo.DistP(pz, p)
    return pz.r

def FindClosestMatchRot(ang):
    theta = math.radians(ang)
    st, ct = math.sin(theta), math.cos(theta)
    rprobepoints = [ P3(p.x*ct+p.y*st, p.y*ct-p.x*st, p.z)+centrepoint  for p in cprobepoints ]
    def fun2(x):
        v = P3(x[0], x[1], x[2])
        tprobepoints = [ p + v  for p in rprobepoints ]
        ds = [ MeasureDistanceToModel(p)  for p in tprobepoints ]
        return sum((d-proberad)**2  for d in ds)
    st = time.time()
    g2 = scipy.optimize.minimize(fun2, [0,0,0], method='Powell')
    return float(g2.fun), time.time()-st, P3(g2.x[0], g2.x[1], g2.x[2])
        
######MAIN

probepoints = MakeProbePoints()
centrepoint = sum(probepoints, P3(0,0,0))*(1/len(probepoints))
cprobepoints = [ p - centrepoint  for p in probepoints ]
FindClosestMatchRot(1)

degs = list(range(0, 360, 2))
random.shuffle(degs)
aopt = [ ]
for deg in degs:
    aopt.append((deg, FindClosestMatchRot(deg)))
    print(aopt[-1])
    aopt.sort()
    sendactivity("clearallcontours")
    sendactivity(contours=[[(d*0.1, fcm[0]*0.1, 50)  for d, fcm in aopt ]])


# https://youtu.be/nZXjUqLMgxM?t=26m23s
