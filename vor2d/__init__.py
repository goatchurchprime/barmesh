
# the outer contour geometry for which we build the voronoi diagram
from .vbdrivegeometry        import VBdrivegeometry

# the process class for building the diagram
from .vbconstructionbarmesh  import VBconstructionbarmesh
from .vbconstructiongeometry import VBconstructiongeometry


# the process class for extracting the diagram from the build object into a simple barmesh
from .vbdiagramcreation      import VBdiagramcreation
from .vbdiagramcreation      import InsertColPointNodes

from .vbdiagramfeatures      import PlotVBdiagram, PlotVBpath
from .vbdiagramfeatures      import GetBarComponents, GetBarSequences, BarsequenceToGseq, BarsequenceToGseqChamf, MergeVBplots
from .vbdiagramfeatures      import GetBarCutR

def DFullBuildVB(conts):
    from basicgeo import Partition1, I1
    urg = I1.AbsorbList(p.u  for cont in conts  for p in cont)
    vrg = I1.AbsorbList(p.v  for cont in conts  for p in cont)
    upart = Partition1(urg.lo-5.01, urg.hi+5, 23+1)
    vpart = Partition1(vrg.lo-5.03, vrg.hi+5, 23+3)
    bm = VBconstructionbarmesh()
    bm.BuildRectBarMesh(upart, vpart)

    vbs = VBdrivegeometry(I1(urg.lo-50, urg.hi+50), I1(vrg.lo-50, vrg.hi+50))
    for cont in conts:  vbs.AddPoly(cont)
    vbs.AddPolyDone(I1(upart.lo-1, upart.hi+1), I1(vpart.lo-1, vpart.hi+1))
    vbc = VBconstructiongeometry(vbs, bm)
    vbc.goconstruct(maxDnextvbarstack=160000)
    print("Construct time %f vbscloset %f%% calls=%d" % (vbc.vbstotaltime, vbc.vbsclosesttime*100/vbc.vbstotaltime, vbc.vbsclosestcalls))

    vbd = VBdiagramcreation(vbs, vbc)
    vbd.BuildVBDC()
    vbm = vbd.vbm
    #InsertColPointNodes(vbs, vbm)
    return vbs, vbm
