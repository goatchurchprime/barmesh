from basicgeo import P2, Along, I1
from math import sqrt
from .vbsimplegeometrycalculations import AdjacentBisectCut, RightSideBisectCorner, LineLineBisectCut, LineArcBisectCut, LineLineArcMid, AdjLineArcLine, AdjLineArcArc, LineArcArcMid, LineLineLineMid, ArcArcArcMid, ArcArcBisectCut
from .vbdrivegeometry import RightSideHorizonSlice

def FindBisectVectorNorm(vbdg, bar, nodeback):
    nodefore = bar.GetNodeFore(bar.nodeback == nodeback)
    assert nodeback.izone != nodefore.izone
    assert bar.nodemid is not None
    
    isegfore = vbdg.GetIseg(nodefore.izone)
    vbfore = vbdg.GetVBiseg(isegfore)
    isegback = vbdg.GetIseg(nodeback.izone)
    vbback = vbdg.GetVBiseg(isegback)
    pv = nodefore.p - nodeback.p
    
    assert vbback.DIsWithinProjectionZone(bar.nodemid.p, vbdg.GetVBisegnext(isegback), 1e-5)
    assert vbfore.DIsWithinProjectionZone(bar.nodemid.p, vbdg.GetVBisegnext(isegfore), 1e-5)
    Dvbbackdist = vbback.Dist(bar.nodemid.p)
    Dvbforedist = vbfore.Dist(bar.nodemid.p)
    assert abs(abs(Dvbbackdist) - abs(Dvbforedist)) < 1e-6
    
    aja = vbdg.AreSegsAdjI(isegback, isegfore)
    if aja != 0:
        assert vbfore.pfrom == vbback.pto or vbback.pfrom == vbfore.pto
        assert vbfore.IsZeroArc() != vbback.IsZeroArc()
        pbdir = vbfore.vperpnorm if vbfore.IsLine() else vbback.vperpnorm
        
    elif vbback.IsLine() and vbfore.IsLine():
        assert abs(abs(Dvbbackdist) - abs(Dvbforedist)) < 1e-6
        # dot of these two doesn't mean anything (obtuse or acute, who cares?)
        #assert (P2.Dot(vbback.vperpnorm, vbfore.vperpnorm) >= 0.0) == (vbbackdist <= 0.0), ("check", vbbackdist, bar.nodemid.p, P2.Dot(vbback.vperpnorm, vbfore.vperpnorm), [(vbback.pfrom, vbback.pto), (vbfore.pfrom, vbfore.pto)])
        if (Dvbbackdist > 0) == (Dvbforedist > 0):
            pbdir = P2.ZNorm((vbback.vperpnorm + vbfore.vperpnorm)*0.5)
        else:
            pbdir = P2.ZNorm((vbback.vperpnorm - vbfore.vperpnorm)*0.5)

    elif vbback.IsLine() or vbfore.IsLine():
        vbline = vbback if vbback.IsLine() else vbfore
        vbarc = vbfore if vbback.IsLine() else vbback
        assert vbarc.IsZeroArc()
        if vbback == vbarc:  # and vice versa
            assert nodeback.r <= abs(vbline.Dist(nodeback.p))+1e-8, (nodeback.r, abs(vbline.Dist(nodeback.p)))
        
        lvecin, lvperpnormdot = (vbline.vperpnorm, vbline.vperpnormdot) if P2.Dot(bar.nodemid.p, vbline.vperpnorm) - vbline.vperpnormdot > 0 else (-vbline.vperpnorm, -vbline.vperpnormdot)
        larcin = bar.nodemid.p - vbarc.pcen
        assert abs(larcin.Len() - bar.nodemid.r) < 1e-6
        assert P2.Dot(bar.nodemid.p, lvecin) - lvperpnormdot > 0, (P2.Dot(bar.nodemid.p, lvecin) - lvperpnormdot)
        assert P2.Dot(vbarc.pcen, lvecin) - lvperpnormdot > 0
        pbdir = P2.ZNorm(P2.CPerp(lvecin) + P2.APerp(larcin)*(1/bar.nodemid.r))

    else:
        assert vbback.IsZeroArc() and vbfore.IsZeroArc()
        pbdir = P2.ZNorm(P2.CPerp(vbfore.pcen - vbback.pcen))
        
    if not pbdir:
        return pbdir
    if P2.Dot(pbdir, P2.CPerp(pv)) < 0.0:
        pbdir = -pbdir
    assert abs(pbdir.Len() - 1.0) < 1e-6
    return pbdir



def FindBisectPoint(vbdg, nodeback, nodefore):
    vbdg.Dfindbisectpoint = (nodeback, nodefore)
    assert nodeback.izone != nodefore.izone

    # first work out the vb segment for the poly the izone corresponds to, decoding contour and element in contour
    isegfore = vbdg.GetIseg(nodefore.izone)
    vbfore = vbdg.GetVBiseg(isegfore)
    isegback = vbdg.GetIseg(nodeback.izone)
    vbback = vbdg.GetVBiseg(isegback)
    
    aja = vbdg.AreSegsAdjI(isegback, isegfore)
    if aja == 1:
        assert vbfore.pfrom == vbback.pto
        assert vbfore.IsZeroArc() != vbback.IsZeroArc()
        bclockwise = vbfore.bclockwise if vbfore.IsZeroArc() else vbback.bclockwise
        pb = AdjacentBisectCut(nodeback.p, nodefore.p, vbfore.prevbisin, vbfore.prevbisout, vbfore.pfrom)
        bAdjSide = RightSideBisectCorner(nodeback.p, nodefore.p, vbfore.pfrom, (vbfore.prevbisout if bclockwise else vbfore.prevbisin))

    elif aja == -1:
        assert vbback.pfrom == vbfore.pto
        assert vbfore.IsZeroArc() != vbback.IsZeroArc()
        bclockwise = vbfore.bclockwise if vbfore.IsZeroArc() else vbback.bclockwise
        pb = AdjacentBisectCut(nodefore.p, nodeback.p, vbback.prevbisin, vbback.prevbisout, vbback.pfrom)
        bAdjSide = RightSideBisectCorner(nodefore.p, nodeback.p, vbback.pfrom, (vbback.prevbisout if bclockwise else vbback.prevbisin))

    elif vbback.IsLine() and vbfore.IsLine():
        pb = LineLineBisectCut(nodeback.p, vbback, nodefore.p, vbfore)
        #assert P2.Dot(vbback.vperpnorm, vbfore.vperpnorm) <= 0.0, "check"

    elif vbback.IsLine() and vbfore.IsArc():
        assert abs((P2.Dot(nodefore.p, vbback.vperpnorm) - vbback.vperpnormdot) - vbback.Dist(nodefore.p)) < 1e-10  # just checking for typos
        if nodefore.r <= abs(P2.Dot(nodefore.p, vbback.vperpnorm) - vbback.vperpnormdot):
            pb = LineArcBisectCut(nodeback.p, vbback, nodefore.p, vbfore)
        else:
            pb = None

    elif vbback.IsArc() and vbfore.IsLine():
        assert abs((P2.Dot(nodeback.p, vbfore.vperpnorm) - vbfore.vperpnormdot) - vbfore.Dist(nodeback.p)) < 1e-10  # just checking for typos
        if nodeback.r <= abs(P2.Dot(nodeback.p, vbfore.vperpnorm) - vbfore.vperpnormdot):
            pb = LineArcBisectCut(nodefore.p, vbfore, nodeback.p, vbback)
        else:
            pb = None

    else:
        pb = ArcArcBisectCut(nodeback.p, vbback, nodefore.p, vbfore)
        
    if not pb:
        return pb, False

    # rshback != rshbacknext means we are within the projection zone of the vbback segment
    rshback = RightSideHorizonSlice(pb, vbback.prevbisin, vbback.prevbisout)
    vbbacknext = vbdg.GetVBisegnext(isegback)
    rshbacknext = RightSideHorizonSlice(pb, vbbacknext.prevbisin, vbbacknext.prevbisout)
    
    rshfore = RightSideHorizonSlice(pb, vbfore.prevbisin, vbfore.prevbisout)
    vbforenext = vbdg.GetVBisegnext(isegfore)
    rshforenext = RightSideHorizonSlice(pb, vbforenext.prevbisin, vbforenext.prevbisout)

    vbback.DIsWithinProjectionZone(pb, vbbacknext, 0.0)  # just testing
    vbfore.DIsWithinProjectionZone(pb, vbforenext, 0.0)
    
    if aja != 0:
        #assert rshback != rshbacknext or rshfore != rshforenext, ("adjfail", pb, ndbsegi, ndfsegi)
        if bAdjSide:
            if not ((rshback != rshbacknext) or (rshfore != rshforenext)):
                vbdg.sendactivity(points=(pb, nodeback.p, vbback.prevbisin, vbfore.prevbisin))
                vbdg.sendactivity(contours=[[nodeback.p, nodefore.p], [vbfore.prevbisin, vbfore.prevbisout], [vbback.prevbisin, vbback.prevbisout]])
                print("bAdjSide", bclockwise, nodeback.p, nodefore.p, vbfore.pfrom, (vbfore.prevbisout if bclockwise else vbfore.prevbisin))
                print((rshback, rshbacknext), (rshfore, rshforenext))
            assert (rshback != rshbacknext) or (rshfore != rshforenext)
            biscontained = True
        else:
            #if bAdjSide:
            #    vbdg.sendactivity(points=(pb, P2.ConvertLZ(nodeback.p), vbback.prevbisin, vbfore.prevbisin))
            #    vbdg.sendactivity(contours=[[P2.ConvertLZ(nodeback.p), P2.ConvertLZ(nodefore.p)], [vbfore.prevbisin, vbfore.prevbisout], [vbback.prevbisin, vbback.prevbisout]])
            #assert not bAdjSide, (bAdjType, pb, P2.ConvertLZ(nodeback.p), P2.ConvertLZ(nodefore.p))
            biscontained = False # round the back
    else:
        biscontained = (rshback != rshbacknext) and (rshfore != rshforenext)
        
    return pb, biscontained

# this tells us about being absolutely certain about what cell the endpoints of the contours belong to
# such segmentation can be achieved as the cells are subdivided down the tree.  assign to cell object colour so we can evaluate down the tree when necessary
def DetectMissingThreeway(vbdg, vbbarcell):
    barnodecuts = vbbarcell.barnodecuts
    izones = [ node.izone  for bar, node in barnodecuts ]
    isegs = [ vbdg.GetIseg(izone)  for izone in izones ]
    vbsegs = [ vbdg.GetVBiseg(iseg)  for iseg in isegs ]
    for izone, vbseg in zip(izones, vbsegs):
        if vbseg.IsLine():
            if vbbarcell.DInteriorness(vbseg.pto) > 0:
                if vbdg.GetIzoneAdvance(izone, 1) not in izones:  # is it missing from our list of zones for this cell
                    return True
    return False
    

# returns pt and r.  Will be -r-1 if non-speculative and outside the range
def MakeCenpt3(vbdg, barnodecuts, bspeccell):
    assert len(barnodecuts) == 3
    izones = [ node.izone  for bar, node in barnodecuts ]
    return MakeCenpt3Z(vbdg, izones, bspeccell, barnodecuts)

def MakeCenpt3Z(vbdg, izones, bspeccell, barnodecutsDir=None):
    assert len(izones) == 3
    assert izones[0] != izones[1] and izones[1] != izones[2] and izones[2] != izones[0]
    isegs = [ vbdg.GetIseg(izone)  for izone in izones ]
    vbsegs = [ vbdg.GetVBiseg(iseg)  for iseg in isegs ] 
    nlines = sum(1  for vbseg in vbsegs  if vbseg.IsLine())
    vbdg.Dvbsegs, vbdg.Disegs, vbdg.DbarnodecutsDir = vbsegs, isegs, barnodecutsDir  # record input for rerunning

    # all these functions should return r as well as cenpt
    bwithinprojzone = True
    if nlines == 3:
        cenpt = LineLineLineMid(vbsegs[0], vbsegs[1], vbsegs[2], barnodecutsDir)
        bwithinprojzone = vbsegs[0].DIsWithinProjectionZone(cenpt, vbdg.GetVBisegnext(isegs[0]), 0) and \
                          vbsegs[1].DIsWithinProjectionZone(cenpt, vbdg.GetVBisegnext(isegs[1]), 0) and \
                          vbsegs[2].DIsWithinProjectionZone(cenpt, vbdg.GetVBisegnext(isegs[2]), 0)

    elif nlines == 2:
        if vbsegs[0].IsArc():
            vbsegs[0], vbsegs[1], vbsegs[2] = vbsegs[1], vbsegs[2], vbsegs[0]
            isegs[0], isegs[1], isegs[2] = isegs[1], isegs[2], isegs[0]
        elif vbsegs[1].IsArc():
            vbsegs[0], vbsegs[1], vbsegs[2] = vbsegs[2], vbsegs[0], vbsegs[1]
            isegs[0], isegs[1], isegs[2] = isegs[2], isegs[0], isegs[1]
        assert vbsegs[0].IsLine() and vbsegs[1].IsLine() and vbsegs[2].IsArc()
        
        aja02 = vbdg.AreSegsAdjI(isegs[0], isegs[2])
        aja12 = vbdg.AreSegsAdjI(isegs[1], isegs[2])
        
        if aja02 != 0 and aja12 != 0 and vbsegs[2].IsZeroArc():
            cenpt = vbsegs[2].pcen
            return cenpt, 0.0  # very special exact zero case
        elif aja02 != 0:
            cenpt = AdjLineArcLine(vbsegs[0], (aja02 == -1), vbsegs[2], vbsegs[1], bspeccell)
        elif aja12 != 0:
            cenpt = AdjLineArcLine(vbsegs[1], (aja12 == -1), vbsegs[2], vbsegs[0], bspeccell)
        else:
            cenpt = LineLineArcMid(vbsegs[0], vbsegs[1], vbsegs[2], bspeccell, barnodecutsDir)
            if cenpt:
                bwithinprojzone = vbsegs[0].DIsWithinProjectionZone(cenpt, vbdg.GetVBisegnext(isegs[0]), 0) and \
                                  vbsegs[1].DIsWithinProjectionZone(cenpt, vbdg.GetVBisegnext(isegs[1]), 0)

    elif nlines == 1:
        if vbsegs[1].IsLine():
            assert vbsegs[2] != vbsegs[0]
            vbsegs[0], vbsegs[1], vbsegs[2] = vbsegs[1], vbsegs[2], vbsegs[0]
            isegs[0], isegs[1], isegs[2] = isegs[1], isegs[2], isegs[0]
        elif vbsegs[2].IsLine():
            vbsegs[0], vbsegs[1], vbsegs[2] = vbsegs[2], vbsegs[0], vbsegs[1]
            isegs[0], isegs[1], isegs[2] = isegs[2], isegs[0], isegs[1]
            assert vbsegs[1] != vbsegs[2]
        assert vbsegs[0].IsLine() and vbsegs[1].IsArc() and vbsegs[2].IsArc()
        assert vbsegs[1] != vbsegs[2]
        
        aja01 = vbdg.AreSegsAdjI(isegs[0], isegs[1])
        aja02 = vbdg.AreSegsAdjI(isegs[0], isegs[2])

        if aja01 != 0:
            cenpt = AdjLineArcArc(vbsegs[0], (aja01 == -1), vbsegs[1], vbsegs[2], bspeccell)
        elif aja02 != 0:
            cenpt = AdjLineArcArc(vbsegs[0], (aja02 == -1), vbsegs[2], vbsegs[1], bspeccell)
        else:
            cenpt = LineArcArcMid(vbsegs[0], vbsegs[2], vbsegs[1], bspeccell, barnodecutsDir)
    else:
        cenpt = ArcArcArcMid(vbsegs[0], vbsegs[2], vbsegs[1])
        
    assert bspeccell or cenpt is not None
    r = (abs(vbsegs[0].Dist(cenpt)) if cenpt else 1e10)
    if not bspeccell and not bwithinprojzone:
        r = -r - 1  # just a way to communicate back for now
    return cenpt, r

