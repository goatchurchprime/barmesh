import sys, random, math
sys.path.append(r"/home/goatchurch/geom3d/barmesh")
from basicgeo import P2, P3, Partition1, Along, Quat
import barmesh, triangleboxing, mainfunctions
from trianglebarmesh import TriangleBarMesh
import implicitareaballoffset, implicitareacyloffset

proberad = 5
Nprobes = 20

sendactivity("clearalltriangles")
#fname = "/home/goatchurch/geom3d/stlfiles/impellor1.stl"
fname = "/home/goatchurch/geom3d/stlfiles/barobox.stl"
tbm = TriangleBarMesh(fname)
sendactivity("triangles", codetriangles=tbm.GetBarMeshTriangles())

ibo = implicitareaballoffset.ImplicitAreaBallOffset(tbm)

proberad = 5
Nprobes = 20


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

probepoints = MakeProbePoints()



# to be used to estimate the accuracy of the positioning
def MeasureDistanceToModel(p):
    pz = barmesh.PointZone(0, 30, None)
    ibo.DistP(pz, p)
    return pz.r
ds = [ MeasureDistanceToModel(p)  for p in probepoints ]
sum((d-proberad)**2  for d in ds)


###################
import scipy.optimize
# now offset and rotate this set of points to somewhere else transform these points to 

# translate by centre of mass of set of probe points
centrepoint = sum(probepoints, P3(0,0,0))*(1/len(probepoints))
cprobepoints = [ p - centrepoint  for p in probepoints ]

# single translation without rotation
def fun2(x):
    v = P3(x[0], x[1], x[2])
    tprobepoints = [ p + v  for p in cprobepoints ]
    ds = [ MeasureDistanceToModel(p)  for p in tprobepoints ]
    return sum((d-proberad)**2  for d in ds)
print(fun2([centrepoint.x,centrepoint.y,centrepoint.z]), "is zero")
g2 = scipy.optimize.minimize(fun2, [0,0,0], method='Powell', options={"maxiter":500}) 
print(g2)
print(g2.x, "=?", centrepoint)


# additionally a rotation about the Z axis
theta = random.uniform(-3.14, 3.14)*1
st, ct = math.sin(theta), math.cos(theta)
rprobepoints = [ P3(p.x*ct+p.y*st, p.y*ct-p.x*st, p.z)  for p in cprobepoints ]

def fun3(x):
    v = P3(x[1], x[2], x[3])
    st, ct = math.sin(x[0]), math.cos(x[0])
    tprobepoints = [ P3(p.x*ct+p.y*st, p.y*ct-p.x*st, p.z)+v  for p in rprobepoints ]
    ds = [ MeasureDistanceToModel(p)  for p in tprobepoints ]
    return sum((d-proberad)**2  for d in ds)
print(fun3([-theta, centrepoint.x,centrepoint.y,centrepoint.z]), "is zero")

# proves it doesn't work well unless only slightly deviated from the chosen angle
g3 = scipy.optimize.minimize(fun3, [-theta-0.2,0,0,0], method='Powell') 
print(g3)
print(g3.x, "=?", -theta, centrepoint)

# option to solve the optimal translation part separately using the direction vectors
# and then do the rotation part separately

# or even plot the value for different rotations to see the number of minima




#############
cprobepoints = [ p - centrepoint  for p in probepoints ]
def fun3(x):
    v = P3(x[0], x[1], x[2])
    tprobepoints = [ p + v  for p in cprobepoints ]
    ds = [ MeasureDistanceToModel(p)  for p in tprobepoints ]
    return sum((d-proberad)**2  for d in ds)
fun2([centrepoint.x,centrepoint.y,centrepoint.z])
g2 = scipy.optimize.minimize(fun2, [0,0,0], method='Powell', options={"maxiter":500}) 
print(g)
print(g.x, "=?", centrepoint)


# create a random rotation and offset
qa = PointsInRadRange(P3(0,0,0), 0, 0.1, 1)[0]
qa = P3(0,0,0)
qr = Quat(math.sqrt(1-qa.Lensq()), qa.x, qa.y, qa.z)
qrv = qr.VecDots()
def PTrans(p):
    p0 = p - centrepoint
    return P3(P3.Dot(p0, qrv[0]), P3.Dot(p0, qrv[1]), P3.Dot(p0, qrv[2]))
rprobepoints = [ PTrans(p)  for p in probepoints ]
sendactivity(points=rprobepoints, materialnumber=1)


###################
print(sum(rprobepoints, P3(0,0,0)))  # WLOG centre is origin, which we can rotate about

import scipy.optimize
ds = [ MeasureDistanceToModel(p)  for p in rprobepoints ]
sum((d-proberad)**2  for d in ds)

def fun(x):
    q = Quat(x[0], x[1], x[2], x[3])
    v = P3(x[4], x[5], x[6])
    qv = q.VecDots()
    tprobepoints = [ P3(P3.Dot(p, qv[0]), P3.Dot(p, qv[1]), P3.Dot(p, qv[2])) + v  for p in rprobepoints ]
    ds = [ MeasureDistanceToModel(p)  for p in tprobepoints ]
    return sum((d-proberad)**2  for d in ds)


g = scipy.optimize.minimize(fun, [1,0,0,0,0,0,0], method='Powell') 
print(g)

fun([1,0,0,0,centrepoint.x,centrepoint.y,centrepoint.z])


