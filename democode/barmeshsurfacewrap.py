import sys
from basicgeo import P2, P3, Partition1, Along
import barmesh, implicitareaballoffset, implicitareacyloffset
import mainfunctions
from barmeshslicer import BarMeshSlicer
from tribarmes import TriangleBarMesh, MakeTriangleBoxing
from barmeshinits import BuildCubeBarMesh, GetAllBarNodeFaces, StarSplitFace

sendactivity("clearalltriangles")        
sendactivity("clearallcontours")
fname = "/home/goatchurch/geom3d/barmesh/stlsamples/frameguide.stl"
tbm = TriangleBarMesh(fname)
ibo = implicitareaballoffset.ImplicitAreaBallOffset(tbm)
sendactivity("triangles", codetriangles=tbm.GetBarMeshTriangles())

bm = barmesh.BarMesh()
xlo, xhi, ylo, yhi, zlo, zhi = -1, 1, -2, 1, 0, 3
rad = 2.5
rex = 7
xlo, xhi, ylo, yhi, zlo, zhi = tbm.xlo-rex, tbm.xhi+rex, tbm.ylo-rex, tbm.yhi+rex, tbm.zlo-rex, tbm.zhi+rex
BuildCubeBarMesh(bm, xlo, xhi, ylo, yhi, zlo, zhi)
#bm.PlotCellLinks(sendactivity)

for bar, node in GetAllBarNodeFaces(bm):
    StarSplitFace(bm, bar, node)
sendactivity(contours=[(bar.nodeback.p, bar.nodefore.p)  for bar in bm.bars  if not bar.bbardeleted])
lll

bm.PlotCellLinks(sendactivity)

#for node in bm.nodes:
#    node.pointzone.v = P3.ZNorm(P3(node.pointzone.v[0], node.pointzone.v[1], node.pointzone.v[2]*1.5))
for node in bm.nodes:
    node.pointzone.bsurfacecontact = False

def ResetVectorsInwards():
    nodebarmap = { }
    for bar in bm.bars:
        if not bar.bbardeleted:
            nodebarmap.setdefault(bar.nodeback, bar)
            nodebarmap.setdefault(bar.nodefore, bar)
    assert len(nodebarmap) == len(bm.nodes)
    for node, bar in nodebarmap.items():
        lbar = bar
        surrnodes = [ ]
        Dcounter = 0
        while True:
            surrnodes.append(lbar.GetNodeFore(lbar.nodeback == node))
            lbar = lbar.GetForeRightBL(lbar.nodefore == node)
            if lbar == bar:
                break
            Dcounter += 1
            assert Dcounter < 1000, "Infinite loop in GetBarForeLeft"
        if len(surrnodes) == 2:
            continue
        vs = [surrnode.p - node.p  for surrnode in surrnodes]
        cvs = [P3.Cross(vs[i+1  if i != len(vs)-1  else 0], vs[i])  for i in range(len(vs))]
        newv = P3.ZNorm(sum(cvs, P3(0,0,0)))
        node.pointzone.v = newv
    
        
def PrepProjectIn(d):
    nprojs = 0
    for node in bm.nodes:
        if not node.pointzone.bsurfacecontact:
            vp = node.pointzone.v * d
            lam = ibo.Cutpos(node.p, vp, None, rad)
            node.newp = node.p + vp*min(1,lam)
            node.newbsurfacecontact = (lam < 1)
            nprojs += 1
    return nprojs

def BarCollides():
    for bar in bm.bars:
        bar.lamcollide = None
        if bar.bbardeleted:  continue
        if bar.nodeback.pointzone.bsurfacecontact or bar.nodefore.pointzone.bsurfacecontact:
            continue
        p0 = bar.nodeback.newp if not bar.nodeback.newbsurfacecontact else bar.nodeback.p
        p1 = bar.nodefore.newp if not bar.nodefore.newbsurfacecontact else bar.nodefore.p
        lam0 = ibo.Cutpos(p0, p1-p0, None, rad)
        if lam0 < 1:
            lam1 = ibo.Cutpos(p1, p0-p1, None, rad)
            assert lam1 < 1
            bar.lamcollide = (lam0 + (1-lam1))*0.5
        
        
def CopyProjectIn():
    for node in bm.nodes:
        if not node.pointzone.bsurfacecontact:
            node.p = node.newp
            node.pointzone.bsurfacecontact = node.newbsurfacecontact
            
#bm.PlotCellLinks(sendactivity)

sendactivity(contours=[(bar.nodeback.p, bar.nodefore.p)  for bar in bm.bars  if not bar.bbardeleted])

# the loop section
ResetVectorsInwards()
PrepProjectIn(1.0)
BarCollides()
sendactivity(points=[Along(bar.lamcollide, bar.nodeback.p, bar.nodefore.p)  for bar in bm.bars  if bar.lamcollide is not None])
CopyProjectIn()
sendactivity(contours=[(bar.nodeback.p, bar.nodefore.p)  for bar in bm.bars  if not bar.bbardeleted])

        
def BarProjectInterpolate():
    conts = [ ]
    for bar in bm.bars:
        cont = [ bar.nodeback.p ]
        for i in range(1, 20):
            lam = i/20
            p = Along(lam, bar.nodeback.p, bar.nodefore.p)
            v = P3.ZNorm(Along(lam, bar.nodeback.pointzone.v, bar.nodefore.pointzone.v))
            vp = v * 40
            lam = ibo.Cutpos(p - vp, vp*2, None, rad)
            pp = p + vp*((min(1,lam)-0.5)*2)
            cont.append(pp)
        cont.append(bar.nodefore.p)
        conts.append(cont)
    sendactivity(contours=conts)


# we're going to try to keep the mesh triangular and avoid the model impinging on its edges
# or its faces as it sinks inwards from free space.
# function to look for closest approach on a triangle, not just a point and an edge 
# so that it eventually settles at a distance r+rex from the model
# and can project 1:1 inwards from that value, filling in the whole machining area.  
# the projection vector is from the triangle corners so it doesn't have gaps.  

# there are connected sequences of bars that are in contact with the surface forming a closed contour that traps a pocket
# the sequence building out the mesh forms a retraction sequence defining the volume fill, perhaps
# we will also have the sequence adaptively following and clearing this surface and it's offset down
# but with the steerage done based on the actual cutter locations projected and the projection surface, 
# rather than this drive surface.  

# then start projecting inwards to the object
# then start subdividing where there are discontinuities
# then subdivide areas and keep it up until it is conformal.  
# isolated breakthrough islands where we hit nothing.  these will be joined long term
