from basicgeo import P2, Along, Partition1, I1
from .vbsimplegeometrycalculations import RightSideHorizonSlice, DHorizonSliceDist
import math

sendactivity = None

def ProjToHorizonBox(p, v, urg, vrg):
    assert urg.ContainsStrict(p.u) and vrg.ContainsStrict(p.v)
    assert v.Lensq() != 0.0
    if v.v > 0.0:
        vlam = (vrg.hi - p.v)/v.v
        vpp = P2(p.u + v.u*vlam, vrg.hi)
    elif v.v < 0.0:
        vlam = (vrg.lo - p.v)/v.v
        vpp = P2(p.u + v.u*vlam, vrg.lo)
    if v.u > 0.0:
        ulam = (urg.hi - p.u)/v.u
        upp = P2(urg.hi, p.v + v.v*ulam)
    elif v.u < 0.0:
        ulam = (urg.lo - p.u)/v.u
        upp = P2(urg.lo, p.v + v.v*ulam)
    
    if v.v != 0.0 and (v.u == 0.0 or vlam < ulam):
        pp = vpp
    else:
        pp = upp
    assert urg.Contains(pp.u) and vrg.Contains(pp.v)
    assert (pp.u == urg.lo or pp.u == urg.hi) or (pp.v == vrg.lo or pp.v == vrg.hi) 
    assert abs(P2.Dot(p - pp, P2.APerp(v))) <= 1e-6
    return pp


class VBsegment:
    def __init__(self, pfrom, pto, pcen):
        self.pfrom = pfrom
        self.pto = pto
        self.pcen = pcen  # =None if this is a line
        #self.bclockwise  # corresponds to convex exists when pcen != None
        self.rad = 0.0
        
        #self.vperpnorm
        #self.vperpnormdot

        # normal bisector at pfrom end point extended beyond the extreme bounding box used for RightSideHorizonSlice()
        self.prevbisin = None  # =ProjToHorizonBox(pfrom, vb.vperpnorm, urg, vrg)
        self.prevbisout = None # =ProjToHorizonBox(pfrom, -vb.vperpnorm, urg, vrg)

    def IsLine(self):
        return self.pcen is None
    def IsArc(self):
        return self.pcen is not None
    def IsZeroArc(self):
        return (self.pcen is not None) and (self.pfrom == self.pto)
    def Dist(self, p):
        if self.pcen != None:
            return (p - self.pcen).Len() - self.rad
        return P2.Dot(p, self.vperpnorm) - self.vperpnormdot
        
    def DistE(self, p):
        if self.pcen != None:
            return (p - self.pcen).Len() - self.rad
        v = self.pto - self.pfrom
        lamv = P2.Dot(p - self.pfrom, v)
        if lamv < 0.0:
            return (p - self.pfrom).Len()
        elif lamv > v.Lensq():
            return (p - self.pto).Len()
        return abs(P2.Dot(p, self.vperpnorm) - self.vperpnormdot)

    def DIsWithinProjectionZone(self, p, vbnext, Depsilon):
        #self = vbsegmentslist[ndbsegilist][ndbsegi]
        #vbnext = vbsegmentslist[ndbsegilist][ndbsegi+1 if ndbsegi+1 != len(vbsegmentslist[ndbsegilist]) else 0]
        assert vbnext.pfrom == self.pto, (vbnext.pfrom, self.pto)
        assert self.pcen or P2.Dot(P2.CPerp(self.prevbisout - self.prevbisin), vbnext.prevbisout - vbnext.prevbisin) < 1e-6
        rsh = RightSideHorizonSlice(p, self.prevbisin, self.prevbisout)
        rshnext = RightSideHorizonSlice(p, vbnext.prevbisin, vbnext.prevbisout)
        assert not self.IsLine() or not (not rsh and rshnext)   # line types are parallel and consistently directionalized oriented
        if rsh != rshnext:
            return True
        elif Depsilon == 0.0:
            return False
        rshDist = abs(DHorizonSliceDist(p, self.prevbisin, self.prevbisout))
        rshnextDist = abs(DHorizonSliceDist(p, vbnext.prevbisin, vbnext.prevbisout))
        return min(rshDist, rshnextDist) < Depsilon

        
class VBdrivegeometry:
    def __init__(self, urg, vrg):
        self.urg = urg
        self.vrg = vrg
        self.vbsegmentslistT = [ ]
        #self.vbsegmentslist = [ ]
        
    def MakeVBlinesegment(self, pfrom, pto):
        vb = VBsegment(pfrom, pto, None)
        vb.vperpnorm = P2.APerp(P2.ZNorm(pto - pfrom))  # points outwards(!)
        vb.vperpnormdot = P2.Dot(pfrom, vb.vperpnorm)
        vb.prevbisin = ProjToHorizonBox(pfrom, vb.vperpnorm, self.urg, self.vrg)
        vb.prevbisout = ProjToHorizonBox(pfrom, -vb.vperpnorm, self.urg, self.vrg)
        return vb
        
    def MakeVBheadnode(self, vbprevlinseg):
        vb = VBsegment(vbprevlinseg.pto, vbprevlinseg.pto, vbprevlinseg.pto)
        vb.prevbisin = ProjToHorizonBox(vbprevlinseg.pto, vbprevlinseg.vperpnorm, self.urg, self.vrg)
        vb.prevbisout = ProjToHorizonBox(vbprevlinseg.pto, -vbprevlinseg.vperpnorm, self.urg, self.vrg)
        return vb

    # insert a line segment and tail node (another zero length VBsegment) between
    def AddPoly(self, pts):
        assert pts[0] == pts[-1]
        vbsegments = [ ]
        for i in range(1, len(pts)):
            assert self.urg.ContainsStrict(pts[i].u) and self.vrg.ContainsStrict(pts[i].v)
            if pts[i-1] != pts[i]:
                vbsegments.append(self.MakeVBlinesegment(pts[i-1], pts[i]))
                vbsegments.append(self.MakeVBheadnode(vbsegments[-1]))
            
        # set the spin values on the zero arcs
        for j in range(1, len(vbsegments), 2):
            assert vbsegments[j].pcen is not None
            assert vbsegments[j].pcen == vbsegments[j].pfrom == vbsegments[j].pto
            j1 = j+1 if j+1 < len(vbsegments) else 0
            assert vbsegments[j1].pcen is None
            assert vbsegments[j1].pfrom == vbsegments[j].pcen
            vh = vbsegments[j].prevbisin - vbsegments[j].prevbisout
            vh1 = vbsegments[j1].prevbisin - vbsegments[j1].prevbisout
            vbsegments[j].bclockwise = P2.Dot(P2.CPerp(vh1), vh) > 0.0
        
        self.vbsegmentslistT.append(vbsegments)
        #sendactivity(contours=[(vb.prevbisin,vb.prevbisout)  for vb in vbs.vbsegmentslist[0]], materialnumber=3)

    def AddPolyDone(self, lurg, lvrg):
        self.vbsegmentslist = self.vbsegmentslistT
        del self.vbsegmentslistT
        minseglength = min(min((vb.pto-vb.pfrom).Len()  for vb in vbsegments  if vb.IsLine())  for vbsegments in self.vbsegmentslist)
        print("minseglength", minseglength)
        self.vbdgboxset = VBdrivegeometryboxset(self.vbsegmentslist, lurg, lvrg)

    def VBclosest(self, pz, p):  # pz = PointZone: (izone, r, v)
        if self.vbdgboxset:
            self.vbdgboxset.VBXclosest(pz, p)

            # for checking the non-boxing case
            #lpz = barmesh.PointZone(-1, -1.0, None)
            #self.VBclosestNoboxing(lpz, p)
            #assert (lpz.r, lpz.izone) == (pz.r, pz.izone), ((lpz.r, lpz.izone), (pz.r, pz.izone), p)
        else:
            self.VBclosestNoboxing(pz, p)
        
    def VBclosestNoboxing(self, pz, p):
        for li, vbsegments in enumerate(self.vbsegmentslist):
            vbprev = vbsegments[-1]
            rshprev = RightSideHorizonSlice(p, vbprev.prevbisin, vbprev.prevbisout)
            iprev = len(vbsegments) - 1
            for i, vb in enumerate(vbsegments):
                rsh = RightSideHorizonSlice(p, vb.prevbisin, vb.prevbisout)
                if rsh != rshprev:   # are we between the extended normal line
                    #Dvbswithinslice.append(vbprev)
                    lr = abs(vbprev.Dist(p))
                    if pz.izone == -1 or lr < pz.r:
                        pz.izone = iprev*len(self.vbsegmentslist) + li
                        pz.r = lr  # and v
                vbprev = vb
                rshprev = rsh
                iprev = i
        
    def MakePointZoneVBS(self, nodes):  # to run in parallel across all nodes
        for node in nodes:
            self.VBclosest(node, node.p)
            
    def AreSegsAdjI(self, iseg0, iseg1):
        if iseg0[0] == iseg1[0]:
            vN = len(self.vbsegmentslist[iseg0[0]])
            if iseg1[1] == (iseg0[1]+1 if iseg0[1]+1 != vN else 0):
                return 1
            if iseg0[1] == (iseg1[1]+1 if iseg1[1]+1 != vN else 0):
                return -1
        return 0

    def GetIseg(self, izone):
        return (izone % len(self.vbsegmentslist), izone // len(self.vbsegmentslist))
    def GetVBiseg(self, iseg):
        return self.vbsegmentslist[iseg[0]][iseg[1]]
    def GetVBisegnext(self, iseg):
        return self.vbsegmentslist[iseg[0]][iseg[1]+1  if iseg[1]+1 != len(self.vbsegmentslist[iseg[0]])  else 0]

    def GetIzoneAdvance(self, izone, adv):
        assert adv == -1 or adv == 1
        iseg0, iseg1 = (izone % len(self.vbsegmentslist), izone // len(self.vbsegmentslist))
        if adv == 1:
            iseg1n = iseg1+1  if iseg1+1 != len(self.vbsegmentslist[iseg0])  else 0
        else:
            iseg1n = iseg1-1  if iseg1 != 0  else len(self.vbsegmentslist[iseg0])-1
        return iseg1n*len(self.vbsegmentslist) + iseg0

    def PlotContours(self):
        return [ [vbseg.pfrom  for vbseg in vbsegments ]  for vbsegments in self.vbsegmentslist ]
        #return [ [vbseg.pfrom  for vbseg in vbsegments  if not vbseg.IsZeroArc()]  for vbsegments in self.vbsegmentslist ]  # skips one edge

    def VBwinding(self, p):
        res = 0
        for vbsegments in self.vbsegmentslist:
            for vb in vbsegments:
                bfromvd = (vb.pfrom.v < p.v)
                btovd = (vb.pto.v < p.v)
                if bfromvd == btovd:
                    continue
                assert vb.IsLine()
                lam = (p.v - vb.pfrom.v) / (vb.pto.v - vb.pfrom.v)
                pu = Along(lam, vb.pfrom.u, vb.pto.u)
                if pu > p.u:
                    res += 1 if vb.pfrom.v < vb.pto.v else -1
        return res
        
def SRightSideHorizonSlice(curg, cvrg, hizin, hizout):
    b00 = RightSideHorizonSlice(P2(curg.lo, cvrg.lo), hizin, hizout)
    b01 = RightSideHorizonSlice(P2(curg.lo, cvrg.hi), hizin, hizout)
    if b01 != b00:
        return 0
    b10 = RightSideHorizonSlice(P2(curg.hi, cvrg.lo), hizin, hizout)
    if b10 != b00:
        return 0
    b11 = RightSideHorizonSlice(P2(curg.hi, cvrg.hi), hizin, hizout)
    if b11 != b00:
        return 0
    return 1 if b00 else -1
        
class VBdrivegeometryboxset:
    def __init__(self, vbsegmentslist, urg, vrg):
        self.vbsegmentslist = vbsegmentslist
        self.urg = urg
        self.vrg = vrg
        boxwidth = max(5, urg.Leng()/50, vrg.Leng()/50)
        self.upart = Partition1(urg.lo, urg.hi, int(math.ceil(urg.Leng()/boxwidth)))
        self.vpart = Partition1(vrg.lo, vrg.hi, int(math.ceil(vrg.Leng()/boxwidth)))
        self.eps = boxwidth*1e-3
        self.boxlookup = { }
        
    def MakeBoxsetuv(self, iu, iv):
        curg = I1(self.upart.vs[iu]-self.eps, self.upart.vs[iu+1]+self.eps)
        cvrg = I1(self.vpart.vs[iv]-self.eps, self.vpart.vs[iv+1]+self.eps)
        cp = P2(curg.Along(0.5), cvrg.Along(0.5))
        lres = [ ]
        for li, vbsegments in enumerate(self.vbsegmentslist):
            vbprev = vbsegments[-1]
            srshprev = SRightSideHorizonSlice(curg, cvrg, vbprev.prevbisin, vbprev.prevbisout)
            iprev = len(vbsegments) - 1
            for i, vb in enumerate(vbsegments):
                srsh = SRightSideHorizonSlice(curg, cvrg, vb.prevbisin, vb.prevbisout)
                if srsh != srshprev or srsh == 0:
                    lr = abs(vbprev.DistE(cp))  # effective distance, not projective distance
                    #izone = iprev*len(self.vbsegmentslist) + li
                    lres.append((lr, li, iprev))
                vbprev = vb
                srshprev = srsh
                iprev = i
        lres.sort()
        
        # use diagonal, not half diagonal because distance goes up at same speed as other distance goes down
        diagleng = (P2(curg.hi, cvrg.hi) - P2(curg.lo, cvrg.lo)).Len()
        while lres and lres[-1][0] - lres[0][0] > diagleng:
            lres.pop()
        lres.sort(key=lambda X:X[1:])
        return [iprev*len(self.vbsegmentslist) + li  for (lr, li, iprev) in lres]
        
        
    def VBXclosest(self, pz, p):
        iu = self.upart.GetPart(p.u)
        iv = self.vpart.GetPart(p.v)
        if (iu, iv) not in self.boxlookup:
            self.boxlookup[(iu, iv)] = self.MakeBoxsetuv(iu, iv)
        
        licurr = -1
        for izone in self.boxlookup[(iu, iv)]:
            li, iprev = izone % len(self.vbsegmentslist), izone // len(self.vbsegmentslist)
            if li != licurr:
                licurr = li
                vbsegments = self.vbsegmentslist[li]
                iprevcurr = -1
            if iprev != iprevcurr:
                vbprev = vbsegments[iprev]
                rshprev = RightSideHorizonSlice(p, vbprev.prevbisin, vbprev.prevbisout)
            i = iprev+1 if iprev+1 != len(vbsegments) else 0
            vb = vbsegments[i]
            rsh = RightSideHorizonSlice(p, vb.prevbisin, vb.prevbisout)

            if rsh != rshprev:
                lr = abs(vbprev.Dist(p))
                if pz.izone == -1 or lr < pz.r:
                    pz.izone = iprev*len(self.vbsegmentslist) + li
                    pz.r = lr  # and v
            vbprevcurr = vb
            rshprev = rsh
            iprev = i

