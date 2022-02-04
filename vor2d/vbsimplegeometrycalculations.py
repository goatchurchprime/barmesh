from basicgeo import P2, Along, I1
from math import sqrt


def RightSideHorizonSlice(p, hizin, hizout):
    if p.v <= hizin.v and p.v <= hizout.v:
        return hizin.u > hizout.u
    if p.v > hizin.v and p.v > hizout.v:
        return hizout.u > hizin.u
    lam = (p.v - hizin.v) / (hizout.v - hizin.v)
    hizu = Along(lam, hizin.u, hizout.u)
    return (hizin.v > hizout.v) == (p.u >= hizu)

def DHorizonSliceDist(p, hizin, hizout):
    vh = hizin - hizout
    return P2.Dot(P2.CPerp(vh), p - hizout) / vh.Len()

DLineLineBisectCut = None
def LineLineBisectCut(pi, vbedgei, po, vbedgeo):
    global DLineLineBisectCut
    DLineLineBisectCut = (pi, vbedgei, po, vbedgeo)
    assert vbedgei.pcen is None and vbedgeo.pcen is None
    deipi = vbedgei.Dist(pi)
    deopi = vbedgeo.Dist(pi)
    deipo = vbedgei.Dist(po)
    deopo = vbedgeo.Dist(po)
    # assert abs(deipi) <= abs(deopi) and abs(deopo) <= abs(deipo)  # not true if we are outside their zone

    # solve: Along(lam, deipi, deipo) == Along(lam, deopi, deopo)
    #  deipi + (deipo - deipi) lam = deopi + (deopo - deopi) lam
    #  deipi - deopi =  (deopo - deopi - deipo + deipi) lam
    if (deipi > 0) == (deopi > 0):
        lnum = deipi - deopi
        ldem = deopo - deipo + lnum
    # solve: Along(lam, deipi, deipo) == -Along(lam, deopi, deopo)
    #  deipi + (deipo - deipi) lam = -deopi - (deopo - deopi) lam
    #  deipi + deopi =  (-deopo + deopi - deipo + deipi) lam
    else:
        lnum = deipi + deopi
        ldem = -deopo - deipo + lnum
        

    if ldem == 0.0:
        return None
    lam = lnum / ldem
    if not 0.0 < lam < 1.0:
        # check if other orientation would work?
        lnum = deipi + deopi
        ldem = -deopo - deipo + lnum
        if ldem == 0.0:
            return None
        lam = lnum / ldem
        if not 0.0 < lam < 1.0:
            return None
        
    pr = Along(lam, pi, po)
    
    #assert abs(vbedgei.Dist(pr) - vbedgeo.Dist(pr)) < 1e-6
    assert abs(abs(vbedgei.Dist(pr)) - abs(vbedgeo.Dist(pr))) < 1e-6
    return pr

DLineArcBisectCutVars = None
def LineArcBisectCut(p0, vbline0, p1, vbarc1):
    global DLineArcBisectCutVars
    DLineArcBisectCutVars = (p0, vbline0, p1, vbarc1)
    assert vbline0.pcen is None and vbarc1.pcen is not None
    assert vbarc1.rad == 0.0

    pv = p1 - p0
    pvsq = pv.Lensq()
    # Dot(pi + pv * clam - vbedgeo.pcen, pv) = 0
    clam = P2.Dot(pv, vbarc1.pcen - p0) / pvsq
    pclos = p0 + pv * clam
    dclossq = (pclos - vbarc1.pcen).Lensq()

    dlclos = P2.Dot(pclos, vbline0.vperpnorm) - vbline0.vperpnormdot
    dvl = P2.Dot(pv, vbline0.vperpnorm)
    vl = dvl / pv.Len()
    # solve dsq = dlam^2 * pvsq + dclossq = (dlclos + dvl * dlam)^2
    qa = dvl**2 - pvsq
    qb2 = dlclos * dvl
    qc = dlclos**2 - dclossq
    assert qa < 1e-10, (dvl**2, pvsq, qa, qb2, qc)  # since vbline0.vperpnorm.Len() == 1, so dvl <= pv.Len()
        
    qdq = qb2**2 - qa*qc
    if qdq < 0.0:
        return None
    if qa != 0.0:
        qs = sqrt(qdq) / qa
        qm = -qb2 / qa
        q0 = qm + qs  # from line closer to arc closer root accounting for qa<0
    else:
        q0 = -qc/(2*qb2)  # place holder

    if qa > -1e-6:
        Dv = q0**2 * qa + q0 * qb2 * 2 + qc
        q0l = -qc/(2*qb2)  # from line closer to arc closer root accounting for qa<0
        Dvl = q0l**2 * qa + q0l * qb2 * 2 + qc
        if abs(Dvl) < abs(Dv):
            q0 = q0l 
    Dv = q0**2 * qa + q0 * qb2 * 2 + qc
    assert abs(Dv) < 1e-5, (Dv, qa, qb2, qc)
    
    lam = clam + q0
    assert 0<lam<1, (lam, q0, qm - qs, clam + qm - qs)
    
    pc = Along(lam, p0, p1)
    Dd0 = vbline0.Dist(pc)
    Dd1 = vbarc1.Dist(pc)
    assert abs(abs(Dd0) - abs(Dd1)) < 1e-6, (Dd0, Dd1, Dd0 - Dd1)
    #print(lam, (Dd0, dlclos + dvl*q0), (Dd1, sqrt(dclossq + q0**2*pvsq)))
    #dlclos = P2.Dot(pclos, vbline0.vperpnorm) - vbline0.vperpnormdot
    return pc

def ArcArcBisectCut(p0, vbarc0, p1, vbarc1):
    assert vbarc0.pcen is not None and vbarc1.pcen is not None
    assert vbarc0.rad == 0.0 and vbarc1.rad == 0.0  # no real arcs for now
    vcen = vbarc1.pcen - vbarc0.pcen
    vcensq = vcen.Lensq()
    dp0 = P2.Dot(p0 - vbarc0.pcen, vcen)
    dp1 = P2.Dot(p1 - vbarc0.pcen, vcen)
    # solve Along(lam, dp0, dp1) = vcensq/2
    if dp0 == dp1:
        return None
    lam = (vcensq/2 - dp0) / (dp1 - dp0)
    if not 0.0 < lam < 1.0:
        return None
    pc = Along(lam, p0, p1)
    Dd0 = vbarc0.Dist(pc)
    Dd1 = vbarc1.Dist(pc)
    assert abs(Dd0 - Dd1) < 1e-6
    return pc

def LineLineMidLine(vbline0, vbline1, s0, s1):
    assert abs(s0) == 1 and abs(s1) == 1
    #vbline0.vperpnorm . pi = vbline0.vperpnormdot
    #vbline1.vperpnorm . pi = vbline1.vperpnormdot
    a, b, c, d = vbline0.vperpnorm.u*s0, vbline0.vperpnorm.v*s0, vbline1.vperpnorm.u*s1, vbline1.vperpnorm.v*s1
    u, v = vbline0.vperpnormdot*s0, vbline1.vperpnormdot*s1
    det = a*d - b*c
    pi = P2((d*u - b*v)/det, (-c*u + a*v)/det)       # distance 0
    assert abs(P2.Dot(pi, vbline0.vperpnorm) - vbline0.vperpnormdot) <= 1e-6
    assert abs(P2.Dot(pi, vbline1.vperpnorm) - vbline1.vperpnormdot) <= 1e-6
    u1, v1 = u+1, v+1
    pi1 = P2((d*u1 - b*v1)/det, (-c*u1 + a*v1)/det)  # distance 1
    vi1 = pi1 - pi
    return pi, vi1

# solve vbline2.Dist(vbline.pi + vi1*d) == d
def MidLineToLine(pi, vi1, vbline2, s2):
    di0 = P2.Dot(pi, vbline2.vperpnorm)*s2 - vbline2.vperpnormdot*s2
    dv1 = P2.Dot(vi1, vbline2.vperpnorm)*s2
    # Dot(pi + vi1 * d, vbline2.vperpnorm) - vbline2.vperpnormdot = d
    if 1 - dv1 == 0:
        return 1e10
    d = di0 / (1 - dv1)
    return d
    
def LineLineLineMid(vbline0, vbline1, vbline2, barnodecutsDir):
    assert vbline0.IsLine() and vbline1.IsLine() and vbline2.IsLine()
    # order to favour the gentlest line
    dl01 = abs(P2.Dot(vbline0.vperpnorm, vbline1.vperpnorm))
    dl12 = abs(P2.Dot(vbline1.vperpnorm, vbline2.vperpnorm))
    dl20 = abs(P2.Dot(vbline2.vperpnorm, vbline0.vperpnorm))
    if dl12 < min(dl01, dl20):
        vbline0, vbline1, vbline2 = vbline1, vbline2, vbline0
    if dl20 < min(dl01, dl12):
        vbline0, vbline1, vbline2 = vbline2, vbline0, vbline1
    
    pi, vi1 = LineLineMidLine(vbline0, vbline1, 1, 1)
    d = MidLineToLine(pi, vi1, vbline2, 1)
    
    pd = pi + vi1 * d
    assert abs(vbline0.Dist(pd) - d) <= 1e-6
    assert abs(vbline1.Dist(pd) - d) <= 1e-6, (vbline1.Dist(pd), d)
    assert abs(vbline2.Dist(pd) - d) <= 1e-6
    
    # check if within cell (actually 3 sides of it) and try alternative line orientations in cases when the contours are misoriented
    if barnodecutsDir:
        dvas = [ P2.Dot(pd - barnodeDir[1].p, P2.CPerp(barnodeDir[0].GetNodeFore(barnodeDir[0].nodeback == barnodeDir[1]).p - barnodeDir[1].p))  for barnodeDir in barnodecutsDir]
        if min(dvas) < 0:
            for s0, s1, s2 in [(-1,1,1),(1,-1,1),(-1,-1,1),(1,1,-1),(-1,1,-1),(1,-1,-1),(-1,-1,-1)]:
                pi, vi1 = LineLineMidLine(vbline0, vbline1, s0, s1)
                d = MidLineToLine(pi, vi1, vbline2, s2)
                pd = pi + vi1 * d
                dvas = [ P2.Dot(pd - barnodeDir[1].p, P2.CPerp(barnodeDir[0].GetNodeFore(barnodeDir[0].nodeback == barnodeDir[1]).p - barnodeDir[1].p))  for barnodeDir in barnodecutsDir]
                if min(dvas) > 0:
                    break
    
    return pd
    
def AdjLineArcLine(vbline0, bfromend, vbarc1, vbline2, bspeccell):
    assert vbline0.IsLine() and vbarc1.IsArc() and vbline2.IsLine()
    pc = vbline0.pfrom if bfromend else vbline0.pto
    d = MidLineToLine(pc, vbline0.vperpnorm, vbline2, 1)
    pd = pc + vbline0.vperpnorm * d
    if (d < 0) != vbarc1.bclockwise:
        d = MidLineToLine(pc, -vbline0.vperpnorm, vbline2, 1)
        pd = pc - vbline0.vperpnorm * d
        assert bspeccell or abs(-vbline0.Dist(pd) - d) <= 1e-6, (vbline0.Dist(pd), d)
    else:
        assert bspeccell or abs(vbline0.Dist(pd) - d) <= 1e-6, (vbline0.Dist(pd), d)

    #print(d, vbline0.Dist(pd), vbarc1.Dist(pd), vbline2.Dist(pd))
    if not bspeccell:  # don't worry about precision if it's a speculative cell
        assert abs(vbarc1.Dist(pd) - abs(d)) <= 1e-6, (vbarc1.Dist(pd), d)
        assert abs(vbline2.Dist(pd) - d) <= 1e-6
    
    return pd
    
DLineArcArcMid = None
def LineArcArcMid(vbline0, vbarc1, vbarc2, bspeccell, barnodecutsDir):
    global DLineArcArcMid
    DLineArcArcMid = (vbline0, vbarc1, vbarc2, barnodecutsDir)
    assert vbline0.IsLine() and vbarc1.IsArc() and vbarc2.IsArc()
    va = vbarc2.pcen - vbarc1.pcen
    valensq = va.Lensq()
    ah = (vbarc2.pcen + vbarc1.pcen)*0.5
    
    # pd = ah + CPerp(va) * k
    hld = P2.Dot(vbline0.vperpnorm, ah) - vbline0.vperpnormdot
    dpnpva = P2.Dot(vbline0.vperpnorm, P2.APerp(va))
    # vbline0.Dist(pd) = hld + k * dpnpva
    # vbarc.Dist(pd) = sqrt(valensq/4 + valensq*k**2)
    # solve |vbline0.Dist(pd)| - vbarc.Dist(pd) == 0
    #   (hld + k * dpnpva)**2 - (valensq/4 + valensq*k**2) == 0
    qa = dpnpva**2 - valensq
    qb2 = hld * dpnpva
    qc = hld**2 - valensq/4
    # qa approaches zero when line is parallel to the line joining the two centres
    
    if abs(qa) > 1e-20:
        qdq = qb2**2 - qa*qc
        assert qdq >= -1e-6, (pi, qdq)
        qs = sqrt(max(qdq, 0.0)) / qa
        qm = -qb2 / qa
        k = qm + qs
        k1 = qm - qs
        if qb2 != 0.0:
            k2 = -qc/(2*qb2)
            if abs(qa*k2**2 + 2*qb2*k2 + qc) < abs(qa*k**2 + 2*qb2*k + qc):
                k = k2
    else:
        k = -qc/(2*qb2)
        
    assert abs(qa*k**2 + 2*qb2*k + qc) < 1e-6, ("bad quadratic solution", qa, qb2, qc, k, abs(qa*k**2 + 2*qb2*k + qc))
    
    pd = ah + P2.APerp(va) * k
    d = hld + k * dpnpva

    if barnodecutsDir:
        dvas = [ P2.Dot(pd - barnodeDir[1].p, P2.CPerp(barnodeDir[0].GetNodeFore(barnodeDir[0].nodeback == barnodeDir[1]).p - barnodeDir[1].p))  for barnodeDir in barnodecutsDir]
        if min(dvas) < 0:
            #print("choose other root of LineArcArcMid", dvas)
            #sendactivity(points=[pd], materialnumber=1)
            pd = ah + P2.APerp(va) * k1
            #sendactivity(points=[pd], materialnumber=1)
            #sendactivity(contours=[[vbline0.pfrom, vbline0.pto]], materialnumber=1)
            #sendactivity(points=[vbarc1.pcen, vbarc2.pcen], materialnumber=1)
            #sendactivity(contours=[[node.p, bar.GetNodeFore(bar.nodeback == node).p]  for bar, node in barnodecutsDir], materialnumber=1)
            d = hld + k1 * dpnpva
            dva1s = [ P2.Dot(pd - barnodeDir[1].p, P2.CPerp(barnodeDir[0].GetNodeFore(barnodeDir[0].nodeback == barnodeDir[1]).p - barnodeDir[1].p))  for barnodeDir in barnodecutsDir]
            if not bspeccell:
                assert min(dva1s) >= 0, dva1s

    if not bspeccell:
        assert abs(vbline0.Dist(pd) - d) <= 1e-6
        assert abs(vbarc1.Dist(pd) - abs(d)) <= 1e-6
        assert abs(vbarc2.Dist(pd) - abs(d)) <= 1e-6
    return pd

def det111(a, b, c, d, e, f):
    return (a*d - b*c) - (a*f - e*b) + (c*f - e*d)
def ArcArcArcMid(vbarc0, vbarc1, vbarc2):
    c0sq, c1sq, c2sq = vbarc0.pcen.Lensq(), vbarc1.pcen.Lensq(), vbarc2.pcen.Lensq()
    sx = det111(c0sq, vbarc0.pcen.v, c1sq, vbarc1.pcen.v, c2sq, vbarc2.pcen.v)
    sy = det111(vbarc0.pcen.u, c0sq, vbarc1.pcen.u, c1sq, vbarc2.pcen.u, c2sq)
    a = det111(vbarc0.pcen.u, vbarc0.pcen.v, vbarc1.pcen.u, vbarc1.pcen.v, vbarc2.pcen.u, vbarc2.pcen.v)
    c = P2(sx, sy)*(0.5/a)
    Dd0, Dd1, Dd2 = (vbarc0.pcen - c).Len(), (vbarc1.pcen - c).Len(), (vbarc2.pcen - c).Len()
    assert abs(Dd0 - Dd1) < 1e-5 and abs(Dd0 - Dd2) < 1e-5
    return c
    
    
DAdjLineArcArc = None
def AdjLineArcArc(vbline0, bfromend, vbarc1, vbarc2, bspeccell):
    global DAdjLineArcArc
    DAdjLineArcArc = (vbline0, bfromend, vbarc1, vbarc2, bspeccell)
    assert vbline0.IsLine() and vbarc1.IsArc() and vbarc2.IsArc()
    assert vbarc1 != vbarc2
    pc = vbline0.pfrom if bfromend else vbline0.pto
    assert pc == vbarc1.pcen

    vi1 = vbline0.vperpnorm

    # Dot(pc + vi1 * dclos - vbarc2.pcen, vi1) = 0
    # Dot(pc - vbarc2.pcen, vi1) + vi1sq * dclos = 0
    assert abs(vi1.Len() - 1.0) < 1e-6
    dclos = P2.Dot(vbarc2.pcen - pc, vi1)
    if dclos == 0.0:
        assert bspeccell
        return None
    pclos = pc + vi1 * dclos
    pclossq = (pclos - vbarc2.pcen).Lensq()
    
    # (dclos + dd)**2 = pclossq + (dd*vi1.Len())**2
    # dclos**2 + 2*dclos*dd = pclossq
    qb2 = -dclos
    qc = pclossq - dclos**2
    dd = -qc / (2*qb2)
    d = dclos + dd
    
    if abs(d) > 1e10:
        assert bspeccell
        return None

    pd = pc + vbline0.vperpnorm * d
    #print(d, vbline0.Dist(pd), vbarc1.Dist(pd), vbline2.Dist(pd))
    assert abs(vbline0.Dist(pd) - d) <= 1e-6
    assert abs(vbarc1.Dist(pd) - abs(d)) <= 1e-6, (vbarc1.Dist(pd), d)
    assert abs(vbarc2.Dist(pd) - abs(d)) <= 1e-6, (vbarc2.Dist(pd), d)
    return pd
    
    
    

DLineLineArcMid = None
def LineLineArcMid(vbline0, vbline1, vbarc2, bspeccell, barnodecutsDir):
    global DLineLineArcMid
    DLineLineArcMid = (vbline0, vbline1, vbarc2)
    assert vbline0.pcen is None and vbline1.pcen is None and vbarc2.pcen is not None
    assert vbarc2.rad == 0.0
    pi, vi1 = LineLineMidLine(vbline0, vbline1, 1, 1)
    if P2.Dot(vbarc2.pcen - pi, vi1) < 0.0:
        vi1 = -vi1
        
    #d = MidLineToArc(pi, vi1, vbarc2)
    # Dot(pi + vi1 * dclos - vbarc2.pcen, vi1) = 0
    # Dot(pi - vbarc2.pcen, vi1) + vi1sq * dclos = 0
    vi1sq = vi1.Lensq()
    dclos = P2.Dot(vbarc2.pcen - pi, vi1) / vi1sq
    pclos = pi + vi1 * dclos
    pclossq = (pclos - vbarc2.pcen).Lensq()
    
    # (dclos + dd)**2 = pclossq + (dd*vi1.Len())**2
    qa = vi1sq - 1
    qb2 = -dclos
    qc = pclossq - dclos**2
    
    assert qa > 0.0, (qa, qb2, qc)
    assert qa > -1e-10, (qa, qb2, qc)
    qdq = qb2**2 - qa*qc

    if qdq < 0.0 and bspeccell:
        return None

    assert qdq >= -1e-6, (pi, qdq)
    qs = sqrt(max(qdq, 0.0)) / qa
    qm = -qb2 / qa

    dd1 = qm + qs  # other solution
    dd = qm - qs   # first (closest) solution
    assert abs(qa*dd**2 + 2*qb2*dd + qc) < 1e-6
    d = dclos + dd
    p = pi + vi1 * d
    
    if barnodecutsDir:
        dvas = [ P2.Dot(p - barnodeDir[1].p, P2.CPerp(barnodeDir[0].GetNodeFore(barnodeDir[0].nodeback == barnodeDir[1]).p - barnodeDir[1].p))  for barnodeDir in barnodecutsDir]
        if min(dvas) < 0:
            print(p, "choose other root of LineLineArcMid", dvas)
            #sendactivity(points=[p], materialnumber=1)
            #sendactivity(points=[p], materialnumber=1)
            #sendactivity(contours=[[vbline0.pfrom, vbline0.pto], [vbline1.pfrom, vbline1.pto]], materialnumber=1)
            #sendactivity(points=[vbarc2.pcen], materialnumber=1)
            #sendactivity(contours=[[node.p, bar.GetNodeFore(bar.nodeback == node).p]  for bar, node in barnodecutsDir], materialnumber=1)
            dd = dd1
            assert abs(qa*dd**2 + 2*qb2*dd + qc) < 1e-6
            d = dclos + dd
            p = pi + vi1 * d
            dva1s = [ P2.Dot(p - barnodeDir[1].p, P2.CPerp(barnodeDir[0].GetNodeFore(barnodeDir[0].nodeback == barnodeDir[1]).p - barnodeDir[1].p))  for barnodeDir in barnodecutsDir]
            if not bspeccell:
                assert min(dva1s) >= 0, (p, dva1s)
    
    return p

def RightSideBisectCorner(pb, pf, hizside, hcorner):
    v = pf - pb  # could use ProjToHorizonBox() then RightSideHorizonSlice() to determin this unambiguously
    vh = hcorner - pb
    vhs = hizside - pb
    rh = P2.Dot(vh, P2.APerp(v))
    rhs = P2.Dot(vhs, P2.APerp(v))
    return (rh < 0.0) != (rhs < 0.0)

def AdjacentBisectCut(pb, pf, hizin, hizout, Dhcorner):
    assert RightSideHorizonSlice(pf, hizin, hizout)
    assert not RightSideHorizonSlice(pb, hizin, hizout)
    assert abs(DHorizonSliceDist(Dhcorner, hizin, hizout)) < 1e-6

    if abs(hizin.u - hizout.u) > abs(hizin.v - hizout.v):
        dhu = hizout.u - hizin.u
        lamf = (pf.u - hizin.u) / dhu
        pfhv = Along(lamf, hizin.v, hizout.v)
        assert (abs(pf.v - pfhv) < 1e-6) or ((pf.v < pfhv) == (dhu < 0.0))
        lamb = (pb.u - hizin.u) / dhu
        pbhv = Along(lamb, hizin.v, hizout.v)
        assert (abs(pb.v - pbhv) < 1e-6) or ((pb.v > pbhv) == (dhu < 0.0)) 
        df = pf.v - pfhv
        db = pb.v - pbhv
        lam = -db / (df - db)
    else:
        dhv = hizout.v - hizin.v
        lamf = (pf.v - hizin.v) / dhv
        pfhu = Along(lamf, hizin.u, hizout.u)
        assert (abs(pf.u - pfhu) < 1e-6) or ((pf.u > pfhu) == (dhv < 0.0))
        lamb = (pb.v - hizin.v) / dhv
        pbhu = Along(lamb, hizin.u, hizout.u)
        assert (abs(pb.u - pbhu) < 1e-6) or ((pb.u < pbhu) == (dhv < 0.0))
        df = pf.u - pfhu
        db = pb.u - pbhu
        lam = -db / (df - db)
    lam = max(0, min(1, lam))  # I1unit.PushInto(lam)
    return Along(lam, pb, pf)



