from basicgeo import P2, P3, Along
from math import sqrt



def PlotVBpath(vbdrivegeometry, bar, rc, rcmax):
    if bar.cellmarkright is None or bar.cellmarkleft is None:
        return [ ]
    izoneleft = bar.cellmarkleft.cellcolour
    izoneright = bar.cellmarkright.cellcolour
    
    isegleft = vbdrivegeometry.GetIseg(izoneleft)
    isegright = vbdrivegeometry.GetIseg(izoneright)
    
    vbsegleft = vbdrivegeometry.GetVBiseg(isegleft)
    vbsegright = vbdrivegeometry.GetVBiseg(isegright)

    if bar.nodeback.pointzone.r < bar.nodefore.pointzone.r:
        nodeback, nodefore = bar.nodeback, bar.nodefore
    else:
        nodeback, nodefore = bar.nodefore, bar.nodeback

    if rc != -1 and nodeback.pointzone.r < rc:
        p0 = GetBarCutR(vbdrivegeometry, bar, rc)
        r0 = rc
    else:
        p0 = nodeback.p
        r0 = nodeback.pointzone.r
    if rcmax != -1 and p0[2] >= rcmax:
        return [ ]
    p0 = P3.ConvertCZ(p0, r0)
    p1 = P3.ConvertCZ(nodefore.p, nodefore.pointzone.r)
    assert rcmax == -1 or p0[2] <= rcmax, (rc, p0[2], rcmax)
    if rcmax != -1 and p1[2] > rcmax:
        p1 = P3.ConvertCZ(GetBarCutR(vbdrivegeometry, bar, rcmax), rcmax)
        
    # linear cases 
    if vbdrivegeometry.AreSegsAdjI(isegleft, isegright) != 0:
        return [p0, p1]
    if vbsegleft.IsLine() and vbsegright.IsLine():
        return [p0, p1]
    #if vbsegleft.IsZeroArc() and vbsegright.IsZeroArc():
    #    return [p0, p1]
    
    # parabolic shape (in xy or z)
    res = [ p0 ]
    for i in range(1, 10):   # should work to tolerance
        r = Along(i/10, r0, p1[2])
        res.append(P3.ConvertCZ(GetBarCutR(vbdrivegeometry, bar, r), r))
    res.append(p1)
    return res
    

def GetBarCutR(vbdrivegeometry, bar, rc):
    if rc == 0:
        if bar.nodeback.pointzone.r == 0:
            return bar.nodeback.p
        assert bar.nodefore.pointzone.r == 0
        return bar.nodefore.p
    assert (bar.nodeback.pointzone.r < rc) != (bar.nodefore.pointzone.r < rc)
    if bar.cellmarkright is None or bar.cellmarkleft is None:
        return (bar.nodeback.p + bar.nodefore.p)*0.5  # should be set using the type of zone we are in and there needs to be MakeColpointNode split points on these lines as well!
    lam = (rc - bar.nodeback.pointzone.r) / (bar.nodefore.pointzone.r - bar.nodeback.pointzone.r)
    assert 0.0 <= lam <= 1.0
        
    izoneleft = bar.cellmarkleft.cellcolour
    izoneright = bar.cellmarkright.cellcolour
    
    isegleft = vbdrivegeometry.GetIseg(izoneleft)
    isegright = vbdrivegeometry.GetIseg(izoneright)
    
    vbsegleft = vbdrivegeometry.GetVBiseg(isegleft)
    vbsegright = vbdrivegeometry.GetVBiseg(isegright)
    
    if vbdrivegeometry.AreSegsAdjI(isegleft, isegright) != 0:
        return Along(lam, bar.nodeback.p, bar.nodefore.p)
    if vbsegleft.IsLine() and vbsegright.IsLine():
        return Along(lam, bar.nodeback.p, bar.nodefore.p)

    if vbsegleft.IsZeroArc() and vbsegright.IsZeroArc():
        vp = vbsegright.pcen - vbsegleft.pcen
        vplen = vp.Len()
        mp = vbsegleft.pcen + vp*0.5
        vpd = sqrt(rc**2 - (vplen*0.5)**2)
        bpos = P2.Dot(P2.APerp(vp), P2.ConvertLZ(bar.nodefore.p - bar.nodeback.p)) > 0
        bdir = bar.nodefore.pointzone.r > bar.nodeback.pointzone.r
        bfac = 1 if (bpos==bdir) else -1   # should be encoded with the segment in some way which branch
        res = mp + P2.APerp(vp)*(vpd/vplen*bfac)
        assert (vbsegleft.Dist(res) - rc) < 1e-4  # and check is in the zone
        assert (vbsegright.Dist(res) - rc) < 1e-4
        return P3.ConvertGZ(res, rc)
        
    vbline = vbsegleft if vbsegleft.IsLine() else vbsegright
    vbarc = vbsegright if vbsegleft.IsLine() else vbsegleft
    assert vbarc.IsZeroArc()

    avec = P2.CPerp(vbline.vperpnorm)
    parccen = P2(P2.Dot(avec, vbarc.pcen), P2.Dot(vbline.vperpnorm, vbarc.pcen) - vbline.vperpnormdot)
    pback = P2(P2.DotLZ(avec, bar.nodeback.p), P2.DotLZ(vbline.vperpnorm, bar.nodeback.p) - vbline.vperpnormdot)
    pfore = P2(P2.DotLZ(avec, bar.nodefore.p), P2.DotLZ(vbline.vperpnorm, bar.nodefore.p) - vbline.vperpnormdot)
    
    # (parccen.y - y)^2 + (x - parccen.x)^2 = y^2
    # 2*parccen.y*y = parccen.y^2 + (x - parccen.x)^2
    # (parccen.y - r)^2 + (x - parccen.x)^2 = r^2
    xdsq = rc**2 - (abs(parccen.v) - rc)**2
    bpos = P2.Dot(avec, P2.ConvertLZ(bar.nodefore.p - bar.nodeback.p)) > 0
    bdir = bar.nodefore.pointzone.r > bar.nodeback.pointzone.r
    bfac = 1 if (bpos==bdir) else -1   # should be encoded with the segment in some way which branch
    x = sqrt(xdsq)*bfac + parccen.u
    bvdir = -1 if parccen.v < 0 else 1 
    res = avec*x + vbline.vperpnorm*(rc*bvdir + vbline.vperpnormdot)
    return P3.ConvertGZ(res, rc)


def MergeVBplots(conts):
    conts = [tuple(cont)  for cont in conts  if cont]
    cmap = { }
    for cont in conts:
        p0, p1 = cont[0], cont[-1]
        if p0 not in cmap:
            cmap[p0] = [ ]
        cmap[p0].append(cont)
        if p1 not in cmap:
            cmap[p1] = [ ]
        cmap[p1].append(cont)

    for jcont in conts:
        for ee in [-1, 0]:
            assert jcont in cmap[jcont[ee]]
        
    sconts = set(conts)
    jconts = [ ]
    while sconts:
        cont = sconts.pop()
        jcont = list(cont[:])
        for ee in [-1, 0]:
            icont = cont
            while True:
                assert jcont[ee] in [icont[0], icont[-1]]
                cconts = cmap[jcont[ee]]
                if len(cconts) != 2:
                    break
                icont = cconts[1-cconts.index(icont)]
                if icont not in sconts:
                    break
                sconts.remove(icont)
                if (icont[0] == jcont[ee]) == (ee == -1):
                    rcont = icont
                else:
                    rcont = list(reversed(icont))
                if ee == -1:
                    jcont.extend(rcont[1:])
                else:
                    jcont = list(rcont)[:-1]+jcont
        if cont in sconts:
            sconts.remove(cont)
        jconts.append(jcont)
    return jconts
    

def PlotVBdiagram(vbdrivegeometry, vbm, rc, coslimit=1, rcmax=-1):
    if rc != -1:
        barcomponents = GetBarComponents(vbm, rc)
        bars = [ ]
        for barcomponent in barcomponents:
            barmax = max(barcomponent, key=lambda X:max(X.nodeback.pointzone.r, X.nodefore.pointzone.r))
            nodemax = barmax.GetNodeFore(barmax.nodefore.pointzone.r > barmax.nodeback.pointzone.r)
            if vbdrivegeometry.VBwinding(P2.ConvertLZ(nodemax.p)) % 2 != 0:
                bars.extend(bar  for bar in barcomponent  if bar.cellmarkright is not None and bar.cellmarkleft is not None)
    else:
        bars = [ bar  for bar in vbm.bars  if not bar.bbardeleted and bar.cellmarkright is not None and bar.cellmarkleft is not None ]
        
    conts = [ ]
    for bar in bars:
        if rc != -1 and coslimit < 1:
            lrc = GetBarRCOS(vbdrivegeometry, bar, rc, coslimit)
            if lrc == -1:
                continue
        else:
            lrc = rc
        conts.append(PlotVBpath(vbdrivegeometry, bar, lrc, rcmax))
        
        # for plotting in the VDirs
        #vbsegleft = vbdrivegeometry.GetVBiseg(vbdrivegeometry.GetIseg(bar.cellmarkleft.cellcolour))
        #vbsegright = vbdrivegeometry.GetVBiseg(vbdrivegeometry.GetIseg(bar.cellmarkright.cellcolour))
        #mp = GetBarCutR(vbdrivegeometry, bar, (bar.nodeback.pointzone.r+bar.nodefore.pointzone.r)*0.5)
        #mp = P2.ConvertLZ(mp)
        #conts.append((mp+GetVDir(vbsegleft, mp)*0.2, mp, mp+GetVDir(vbsegright, mp)*0.2))

    return conts
    



    
def GetBarsRoundNode(bar, node):
    res = [ ]
    lbar = bar
    Dcount = 0
    while True:
        lbar = lbar.GetForeRightBL(lbar.nodefore == node)
        if lbar == None:
            break
        if lbar == bar:
            return res
        res.append(lbar)
        Dcount += 1
        assert Dcount < 100
    lbar = bar
    while True:
        lbar = lbar.GetForeLeftBR(lbar.nodefore == node)
        assert lbar != bar
        if lbar == None:
            break
        res.append(lbar)
        Dcount += 1
        assert Dcount < 100
    return res
            

def GetBarComponents(vbm, rc):
    barset = set(bar  for bar in vbm.bars  if not bar.bbardeleted  if max(bar.nodeback.pointzone.r, bar.nodefore.pointzone.r) > rc)
    barcomponents = [ ]
    while barset:
        sbar = min(barset, key=lambda X:(X.nodeback.i, X.nodefore.i))
        barset.remove(sbar)
        barnodestack = [ ]
        if sbar.nodeback.pointzone.r > rc:
            barnodestack.append((sbar, sbar.nodeback))
        if sbar.nodefore.pointzone.r > rc:
            barnodestack.append((sbar, sbar.nodefore))
        if not barnodestack:
            continue
        barcomponent = [ sbar ]
        while barnodestack:
            bar, node = barnodestack.pop()
            nbars = GetBarsRoundNode(bar, node)
            assert 1 <= len(nbars) <= 2
            for lbar in nbars:
                if lbar in barset:
                    barcomponent.append(lbar)
                    barset.remove(lbar)
                    lnode = lbar.GetNodeFore(lbar.nodeback == node)
                    if lnode.pointzone.r > rc:
                        barnodestack.append((lbar, lnode))
        barcomponents.append(barcomponent)
    return barcomponents


def GetBarSequences(barcomponent, rc):
    barsequences = [ ]
    barset = set(bar  for bar in barcomponent  if min(bar.nodeback.pointzone.r, bar.nodefore.pointzone.r) <= rc)
    while barset:
        cbar = min(barset, key=lambda X:(X.nodeback.i, X.nodefore.i))
        barset.remove(cbar)
        cnode = cbar.GetNodeFore(cbar.nodeback.pointzone.r <= rc)
        barsequence = [ cbar ]
        assert cnode.pointzone.r > rc
        bar, node = cbar, cnode
        while True:
            if bar.GetForeRightBL(bar.nodeback == node) is None:
                break
            assert bar.GetNodeFore(bar.nodeback == node).pointzone.r <= rc
            assert node.pointzone.r > rc
            lbar, lnode = bar, node
            cellmark = lbar.GetCellMarkRightL(lbar.nodefore == lnode)
            Dcount = 0
            while True:
                lbar = lbar.GetForeRightBL(lbar.nodefore == lnode)
                lnode = lbar.GetNodeFore(lbar.nodeback == lnode)
                assert cellmark == lbar.GetCellMarkRightL(lbar.nodefore == lnode)
                assert lbar != bar
                if lnode.pointzone.r <= rc:
                    break
                Dcount += 1
                assert lbar != bar
                assert Dcount < 1000
            bar, node = lbar, lbar.GetNodeFore(lbar.nodeback == lnode)
            if bar == cbar:
                assert bar not in barset
                break
            barset.remove(bar)
            barsequence.append(bar)
        if bar.GetForeRightBL(bar.nodeback == node) is None:
            assert False, "Need to write the iterate reverse direction code"
        else:
            barsequence.append(barsequence[0])  # loop
        barsequences.append(barsequence)
    return barsequences


def BarsequenceToGseq(vbdrivegeometry, barsequence, rc):
    p0 = GetBarCutR(vbdrivegeometry, barsequence[0], rc)
    gseq = [ (p0[0], p0[1], 1) ]
    for bar in barsequence[1:]:
        p1 = GetBarCutR(vbdrivegeometry, bar, rc)
        cellmark = bar.GetCellMarkRightL(bar.nodefore.pointzone.r <= rc)
        iseg = vbdrivegeometry.GetIseg(cellmark.cellcolour)
        vb = vbdrivegeometry.GetVBiseg(iseg)
        if vb.IsZeroArc():
            gseq.append((p1[0], p1[1], 2, vb.pcen[0], vb.pcen[1]))
            assert abs((P2.ConvertLZ(p1)-vb.pcen).Len() - (P2.ConvertLZ(p0)-vb.pcen).Len()) < 1e-5
        else:
            gseq.append((p1[0], p1[1], 1))
            # prove parallel to the line
        p0 = p1
    return gseq


def RChamfProj(p, vb, rc, rout):
    if vb.IsZeroArc():
        assert abs((p - vb.pcen).Len() - rc) < 1e-4
        pC = p - (p - vb.pcen)*(rout/rc)
    else:
        assert vb.IsLine()
        vbd = P2.Dot(p, vb.vperpnorm) - vb.vperpnormdot
        assert abs(abs(vbd) - rc) < 1e-4, (vbd, rc)
        pC = p + vb.vperpnorm*(-rout if vbd > 0 else rout)
    return pC

def GetBarCutRChamf(vbdrivegeometry, bar, rc, rout):
    assert rout <= rc
    cellmarkL = bar.GetCellMarkRightL(bar.nodefore.pointzone.r <= rc)
    cellmarkR = bar.GetCellMarkRightL(bar.nodefore.pointzone.r > rc)
    isegL = vbdrivegeometry.GetIseg(cellmarkL.cellcolour)
    isegR = vbdrivegeometry.GetIseg(cellmarkR.cellcolour)
    vbL = vbdrivegeometry.GetVBiseg(isegL)
    vbR = vbdrivegeometry.GetVBiseg(isegR)
    p = P2.ConvertLZ(GetBarCutR(vbdrivegeometry, bar, rc))
    pL = RChamfProj(p, vbL, rc, rout)
    pR = RChamfProj(p, vbR, rc, rout)
    
    res = [ ]
    if vbL.IsZeroArc():
        res.append((pL[0], pL[1], 2, vbL.pcen[0], vbL.pcen[1]))
    else:
        res.append((pL[0], pL[1], 1))
    
    if vbdrivegeometry.AreSegsAdjI(isegL, isegR) == 0:
        res.append((pR[0], pR[1], 3, p[0], p[1]))
    else:
        assert (pL - pR).Len() < 1e-4
    return res
    
    
def BarsequenceToGseqChamf(vbdrivegeometry, barsequence, rc, rout):
    assert rout <= rc
    #prevbar2 = barsequence[-1]
    #p0 = GetBarCutR(vbdrivegeometry, prevbar, rc)
    gseq = [ ]
    for bar in barsequence:
        lgseq = GetBarCutRChamf(vbdrivegeometry, bar, rc, rout)
        gseq.extend(lgseq)
        continue
    gseq[0] = (gseq[0][0], gseq[0][1], 1)
        
    return gseq



def GetVDir(vbseg, p):
    if vbseg.IsLine():
        if P2.Dot(vbseg.vperpnorm, p - vbseg.pfrom) > 0:
            return -vbseg.vperpnorm
        return vbseg.vperpnorm
        
    assert vbseg.IsZeroArc()
    return P2.ZNorm(vbseg.pcen - p)
    

def GetBarRCOS(vbdrivegeometry, bar, rc, coslimit):
    assert not (bar.nodeback.pointzone.r < rc) or not (bar.nodefore.pointzone.r < rc)
    assert not (bar.cellmarkright is None or bar.cellmarkleft is None)
    
    izoneleft = bar.cellmarkleft.cellcolour
    izoneright = bar.cellmarkright.cellcolour
    
    isegleft = vbdrivegeometry.GetIseg(izoneleft)
    isegright = vbdrivegeometry.GetIseg(izoneright)
    
    vbsegleft = vbdrivegeometry.GetVBiseg(isegleft)
    vbsegright = vbdrivegeometry.GetVBiseg(isegright)
        
    if bar.nodeback.pointzone.r < bar.nodefore.pointzone.r:
        nodeback, nodefore = bar.nodeback, bar.nodefore
    else:
        nodeback, nodefore = bar.nodefore, bar.nodeback
    p1cos = P2.Dot(GetVDir(vbsegleft, P2.ConvertLZ(nodefore.p)), GetVDir(vbsegright, P2.ConvertLZ(nodefore.p)))

    if vbdrivegeometry.AreSegsAdjI(isegleft, isegright) != 0:
        assert abs(p1cos - 1) < 1e-4, p1cos
        return -1
    if coslimit == 1:
        return rc
    p0cos = P2.Dot(GetVDir(vbsegleft, P2.ConvertLZ(nodeback.p)), GetVDir(vbsegright, P2.ConvertLZ(nodeback.p)))
    if vbsegleft.IsLine() and vbsegright.IsLine():
        #assert nodeback.pointzone.r == 0 or abs(p1cos - p0cos) < 1e-4, (p1cos, p0cos, nodeback.p)
        if p1cos < coslimit:
            return rc
        return -1
    assert p0cos <= p1cos, (p0cos, p1cos)
    if p0cos >= coslimit:
        return -1
    if p1cos <= coslimit:
        return rc
    return rc
    
    # to do this limiting properly for parabolas we'd have to set the upper rc on the bar since the 
    # cos value increases with r and will trip through the coslimit.  It's probably not desirable to 
    # have unnecessary disconnections away from the contour that this would make
    return (nodeback.pointzone.r + nodefore.pointzone.r)*0.5

