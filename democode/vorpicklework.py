from basicgeo import P2, P3, Partition1, I1, Along
import barmesh
from vor2d import VBdrivegeometry, VBconstructiongeometry
from vor2d import VBdiagramcreation
from math import tan, radians
from vor2d import PlotVBdiagram, PlotVBpath, GetBarComponents, GetBarCutR, GetBarSequences, BarsequenceToGseq, BarsequenceToGseqChamf, MergeVBplots
from postprocess import PostProcessor

sendactivity("clearalltriangles")
sendactivity("clearallcontours")
sendactivity("clearallpoints")

import pickle
pfile = "/home/goatchurch/geom3d/leathersatchel/vbsvbm8r.pickle"
(vbs, vbm) = pickle.load(open(pfile, "rb"))
sendactivity(contours=vbs.PlotContours(), materialnumber=2)

cleartoolrad = 2*0.5
chamferthickness = 0.1
taperslope = tan(radians(15))
tapertiprad = 0.3*0.5
chamfercleartoolrad = 0.5
chamferz = -chamfercleartoolrad/taperslope + tapertiprad/taperslope
    
cfile = "/home/goatchurch/geom3d/leathersatchel/leather.ngc"
postprocessor = PostProcessor(cfile, sendactivity)

def Main():
  if 1:  # facing tool
    barsequences = MakeOffsets(0)
    OutputGseqsLayer(Facing(4, GetBounds(barsequences[0], 3)), postprocessor, z=0, z1=3, breorder=False)
  postprocessor.zclear = 1
  if 0:  # clear border to depth because stock too big
    barsequences = MakeOffsets(0)
    OutputGseqsLayer(Boxing(GetBounds(barsequences[0], 3)), postprocessor, z=-1.5, z1=1)
  if 0:  # 2mm endmill
    rsteps = [cleartoolrad + chamferthickness + i*cleartoolrad*0.95  for i in range(9)]
    gseqs = Offsets(rsteps)
    postprocessor.zclear = 1
    OutputGseqsLayer(gseqs, postprocessor, z=-1.5, z1=1)
    OutputGseqsLayer(gseqs, postprocessor, z=-3, z1=1)
  if 0:  # tapered cutter
    gseqs = Offsets([chamfercleartoolrad])
    z = -(chamfercleartoolrad-tapertiprad)/taperslope
    OutputGseqsLayer(gseqs, postprocessor, z=z, z1=1)
  if 0:  # tapered cutter
    vvconts = PlotVBdiagram(vbs, vbm, rc=tapertiprad, coslimit=0.9, rcmax=chamfercleartoolrad)
    vconts = MergeVBplots(vvconts)
    OutputVBGseqs(vconts, postprocessor, 1)
    
  sendactivity(contours=postprocessor.contcuts, materialnumber=1)
  sendactivity(contours=postprocessor.contarcs, materialnumber=0)
  sendactivity(contours=postprocessor.contlinks, materialnumber=2)
  postprocessor.Close()




def MakeOffsets(rc):
    barcomponents = GetBarComponents(vbm, rc)
    barsequences = [ ]
    for barcomponent in barcomponents:
        barmax = max(barcomponent, key=lambda X:max(X.nodeback.pointzone.r, X.nodefore.pointzone.r))
        nodemax = barmax.GetNodeFore(barmax.nodefore.pointzone.r > barmax.nodeback.pointzone.r)
        if vbs.VBwinding(P2.ConvertLZ(nodemax.p)) % 2 != 0:
            lbarsequences = GetBarSequences(barcomponent, rc)
            barsequences.extend(lbarsequences)
    return barsequences


def OutputGseqsLayer(gseqs, postprocessor, z, z1, breorder=True):
    gseqs = gseqs[:]
    pf = (0,0)
    postprocessor.zclear = z1
    while gseqs:
        if breorder:
            gseq = min(gseqs, key=lambda X:(X[0][0]-pf[0])**2 + (X[0][1]-pf[1])**2)
        else:
            gseq = gseqs[0]
        gseqs.remove(gseq)
        postprocessor.OutGseq(gseq, z)
        pf = gseq[-1]

    
def OutputVBGseqs(lconts, postprocessor, z1):
    lconts = lconts[:]
    pf = (0,0)
    while lconts:
        lcont = min(lconts, key=lambda X:(X[-1][0]-pf[0])**2 + (X[-1][1]-pf[1])**2)
        lconts.remove(lcont)
        lcont.reverse()
        lcont = [ (p[0], p[1], -(p[2]-tapertiprad)/taperslope)  for p in lcont ]
        postprocessor.OutLcont(lcont)
        pf = lcont[-1]
    print("vblowest", min(p[2]  for cont in postprocessor.contcuts  for p in cont))


def Offsets(rsteps):
    gseqs = [ ]
    for r in rsteps:
        barsequences = MakeOffsets(r)
        if not barsequences:
            break
        lgseqs = [BarsequenceToGseq(vbs, barsequence, r)  for barsequence in barsequences]
        gseqs = lgseqs + gseqs
    return gseqs

def Chamfering(rc, rcout):
    barsequences = MakeOffsets(rc)
    gseqcs = [BarsequenceToGseqChamf(vbs, barsequence, rc, rcout)  for barsequence in barsequences]
    return gseqcs
                
def GetBounds(barsequence, offs):
    pts = [ GetBarCutR(vbs, bar, 0)  for bar in barsequence ] # outer contour
    x0, x1 = min(p[0]  for p in pts)-offs, max(p[0]  for p in pts)+offs
    y0, y1 = min(p[1]  for p in pts)-offs, max(p[1]  for p in pts)+offs
    return x0, x1, y0, y1

    
def Facing(step, bb):
    (x0, x1, y0, y1) = bb
    gseqs = [ ]
    x, y = x0, y0
    while x < x1:
        gseqs.append([(x, y0, 1), (x, y1, 1)])
        x += step
    gseqs.append([(x1, y0, 1), (x1, y1, 1)])
    return gseqs

def Boxing(bb):
    (x0, x1, y0, y1) = bb
    return [[(x0,y0,1), (x1,y0,1), (x1,y1,1), (x0,y1,1), (x0,y0,1)]]

Main()

"""
conts=PlotVBdiagram(vbs, vbm, rc=-1, coslimit=1)
sendactivity(contours=conts, materialnumber=1)
sendactivity(points=[cont[0]  for cont in conts])
sendactivity(points=[cont[-1]  for cont in conts])
"""
