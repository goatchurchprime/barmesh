from basicgeo import P2, P3
import barmesh


class VBDCbmpath:
    #                               izoneL
    # vbbarcell0, i0, vnode0 ... vsbarnodecuts ... vbbarcell1, i1, vnode1   ==> vbar
    #                               izoneR
    def __init__(self, vbbarcell0, i0, vbdc):
        self.vbbarcell0, self.i0 = vbbarcell0, i0
        if vbbarcell0 is None:
            bar, node = vbdc.extcyclevnodes[i0][:2]
            self.vnode0 = vbdc.extcyclevnodes[i0][2]
        else:
            assert len(vbbarcell0.barnodecuts) == 3
            bar, node = vbbarcell0.barnodecuts[i0]
            self.vnode0 = vbdc.vbbarcellvnodemap[vbbarcell0]
        self.izoneL = bar.GetNodeFore(bar.nodefore == node).izone
        self.izoneR = bar.GetNodeFore(bar.nodeback == node).izone
        self.vsbarnodecuts = [ ]
        Dntrack = 0
        vbbarcell, i = vbbarcell0, i0
        while True:
            vbcl = vbdc.barvbcellmap[bar]
            if vbbarcell is not None:
                if len(vbcl) == 1:
                    assert bar in vbdc.extbarcyclemap, (bar.nodeback.p, bar.nodefore.p)
                    vbbarcell, i = None, vbdc.extbarcyclemap[bar]
                    assert vbdc.extcyclevnodes[i][0] == bar
                    self.vnode1 = vbdc.extcyclevnodes[i][2]
                    self.vsbarnodecuts.append((bar, bar.GetNodeFore(bar.nodeback == node)))
                    break
                assert len(vbcl) == 2, "k"
                j = vbcl.index(vbbarcell)
                vbbarcell = vbcl[1-j]
                
            # entry case from ext
            else:
                assert len(vbcl) == 1
                vbbarcell = vbcl[0]  
            node = bar.GetNodeFore(bar.nodeback == node)
            self.vsbarnodecuts.append((bar, node))
            i = vbbarcell.barnodecuts.index((bar, node))
            if len(vbbarcell.barnodecuts) == 3:
                self.vnode1 = vbdc.vbbarcellvnodemap[vbbarcell]
                break
            assert len(vbbarcell.barnodecuts) == 2
            bar, node = vbbarcell.barnodecuts[1-i]
            assert self.izoneL == bar.GetNodeFore(bar.nodefore == node).izone, "k4"
            assert self.izoneR == bar.GetNodeFore(bar.nodeback == node).izone, "k5"

            Dntrack += 1
            assert Dntrack < 1000
        
        # return endpoints and a hash that sorts pairs together
        sk0 = -1-self.i0  if self.vbbarcell0 is None else self.vbbarcell0.GetCellMark().cellcolour*3 + self.i0
        sk = -1-i  if vbbarcell is None else vbbarcell.GetCellMark().cellcolour*3 + i
        self.vbbarcell1, self.i1 = vbbarcell, i
        self.skC = (min(sk0, sk), sk0)

    def DProveReversedMatch(self, ovbdc):
        assert self.skC[0] == ovbdc.skC[0]
        assert (self.izoneL, self.izoneR) == (ovbdc.izoneR, ovbdc.izoneL)
        assert (self.vbbarcell0, self.i0) == (ovbdc.vbbarcell1, ovbdc.i1)
        assert (self.vbbarcell1, self.i1) == (ovbdc.vbbarcell0, ovbdc.i0)
        assert len(self.vsbarnodecuts) == len(ovbdc.vsbarnodecuts)
        assert (self.vnode0, self.vnode1) == (ovbdc.vnode1, ovbdc.vnode0), ((self.vnode0.i, self.vnode1.i), (ovbdc.vnode1.i, ovbdc.vnode0.i))
        return True
        
# build up the nodes and references to nodes
def NewNodeR(vbm, p, r):
    node = vbm.NewNode(P3.ConvertGZ(p, r))
    node.pointzone = barmesh.PointZone(0, r, None)  # we could lose this with r in the point.z position
    return node


class VBdiagramcreation:
    def __init__(self, vbdrivegeometry, vbconstructiongeometry):
        self.vbdrivegeometry = vbdrivegeometry
        self.vbconstructiongeometry = vbconstructiongeometry
        assert vbconstructiongeometry.vbdrivegeometry == vbdrivegeometry
        extbarcycles = GetBarCycles(vbconstructiongeometry.bm)
        assert len(extbarcycles) == 1
        self.extbarcycle = extbarcycles[0]
    
    # separate function from constructor so we can query results after an assertion
    def BuildVBDC(self):
        # build the cellmarks
        self.izonecellmarkmap = dict((izone, barmesh.CellMark(izone))  for izone in (i*len(self.vbdrivegeometry.vbsegmentslist) + li  for li in range(len(self.vbdrivegeometry.vbsegmentslist))  for i in range(len(self.vbdrivegeometry.vbsegmentslist[li]))))
        
        self.vbm = barmesh.BarMesh()
        
        self.extcyclevnodes = [ (bar, node, NewNodeR(self.vbm, bar.nodemid.p, bar.nodemid.r))  for (bar, node) in self.extbarcycle  if bar.nodemid is not None ]
        self.extbarcyclemap = dict((bar, i)  for i, (bar, node, vnode) in enumerate(self.extcyclevnodes))
        vbbarcells3 = [vbbarcell  for vbbarcell in self.vbconstructiongeometry.vbbarcells  if vbbarcell is not None and len(vbbarcell.barnodecuts) == 3]
        # vbbarcells3.sort(key=lambda X:X.cenptR)  # get these nodes in ascending distance order (which can't be satisfied once we insert colpointnodes)
        self.vbbarcellvnodemap = dict((vbbarcell, NewNodeR(self.vbm, vbbarcell.cenpt, vbbarcell.cenptR))  for vbbarcell in vbbarcells3)
                
        # map from bars to cells either side
        self.barvbcellmap = { }   
        for vbbarcell in self.vbconstructiongeometry.vbbarcells:
            if vbbarcell is not None:
                for bar, node in vbbarcell.barnodecuts:
                    if bar not in self.barvbcellmap:
                        self.barvbcellmap[bar] = [ ]
                    self.barvbcellmap[bar].append(vbbarcell)
                    assert len(self.barvbcellmap[bar]) <= 2

        # track through from everywhere to everywhere
        self.allvtrackpairs = [ ]
        for i in range(len(self.extcyclevnodes)):
            self.allvtrackpairs.append(VBDCbmpath(None, i, self))
        for vbbarcell in vbbarcells3:
            for i in range(3):
                self.allvtrackpairs.append(VBDCbmpath(vbbarcell, i, self))
        assert len(self.vbbarcellvnodemap)*3 == len(self.allvtrackpairs) - len(self.extcyclevnodes)

        # sort and merge the tracks
        self.vbbarcellvbarsmap = { }
        self.allvtrackpairs.sort(key=lambda X:X.skC)
        assert (len(self.allvtrackpairs) % 2) == 0
        for i in range(2, len(self.allvtrackpairs), 2):
            assert self.allvtrackpairs[i-1].skC[0] != self.allvtrackpairs[i].skC[0]

        # pick the tracks according to node order for bar creation
        wvtrackpairs = [ ]
        for i in range(1, len(self.allvtrackpairs), 2):
            assert self.allvtrackpairs[i-1].DProveReversedMatch(self.allvtrackpairs[i])
            assert self.allvtrackpairs[i-1].izoneL in self.izonecellmarkmap
            assert self.allvtrackpairs[i-1].izoneR in self.izonecellmarkmap
            if self.allvtrackpairs[i-1].vnode0.i < self.allvtrackpairs[i-1].vnode1.i:
                wvtrackpairs.append(self.allvtrackpairs[i-1])
            else:
                wvtrackpairs.append(self.allvtrackpairs[i])

        # discover the 3-way ordering at each node on the wvtrackpairs pre-bars
        vbbarcellvbar3map = dict((vbbarcell, [None, None, None])  for vbbarcell in vbbarcells3)
        for wvtrackpair in wvtrackpairs:
            if wvtrackpair.vbbarcell0 is not None:
                assert vbbarcellvbar3map[wvtrackpair.vbbarcell0][wvtrackpair.i0] is None
                vbbarcellvbar3map[wvtrackpair.vbbarcell0][wvtrackpair.i0] = wvtrackpair
            if wvtrackpair.vbbarcell1 is not None:
                assert vbbarcellvbar3map[wvtrackpair.vbbarcell1][wvtrackpair.i1] is None
                vbbarcellvbar3map[wvtrackpair.vbbarcell1][wvtrackpair.i1] = wvtrackpair
                
        # construct the internal bars
        for wvtrackpair in wvtrackpairs:
            wvtrackpair.vbar = barmesh.Bar(wvtrackpair.vnode0, wvtrackpair.vnode1)
            wvtrackpair.vbar.cellmarkleft = self.izonecellmarkmap[wvtrackpair.izoneL]
            wvtrackpair.vbar.cellmarkright = self.izonecellmarkmap[wvtrackpair.izoneR]
            
        # construct the external bars with last one reversed to maintain index order
        extvbars = [ barmesh.Bar(self.extcyclevnodes[0][2], self.extcyclevnodes[-1][2]) ] + [ barmesh.Bar(self.extcyclevnodes[i-1][2], self.extcyclevnodes[i][2])  for i in range(1, len(self.extcyclevnodes)) ]

        # construct the links between internal bars and to and from external bars
        for wvtrackpair in wvtrackpairs:
            if wvtrackpair.vbbarcell0 is not None:
                i0b = wvtrackpair.i0 - 1 if wvtrackpair.i0 != 0 else 2
                wvtrackpair.vbar.barbackleft = vbbarcellvbar3map[wvtrackpair.vbbarcell0][i0b].vbar
            else:
                wvtrackpair.vbar.barbackleft = extvbars[wvtrackpair.i0]
                i0ep = wvtrackpair.i0+1 if wvtrackpair.i0+1 != len(self.extcyclevnodes) else 0
                extvbars[i0ep].SetForeRightBL(i0ep == 0, wvtrackpair.vbar)
                if wvtrackpair.i0 == 0:
                    assert extvbars[wvtrackpair.i0].cellmarkright is None
                    extvbars[wvtrackpair.i0].cellmarkright = self.izonecellmarkmap[wvtrackpair.izoneL]
                else:
                    assert extvbars[wvtrackpair.i0].cellmarkleft is None
                    extvbars[wvtrackpair.i0].cellmarkleft = self.izonecellmarkmap[wvtrackpair.izoneL]
                
            if wvtrackpair.vbbarcell1 is not None:
                i1b = wvtrackpair.i1 - 1 if wvtrackpair.i1 != 0 else 2
                wvtrackpair.vbar.barforeright = vbbarcellvbar3map[wvtrackpair.vbbarcell1][i1b].vbar
            else:
                wvtrackpair.vbar.barforeright = extvbars[wvtrackpair.i1]
                i1ep = wvtrackpair.i1+1 if wvtrackpair.i1+1 != len(self.extcyclevnodes) else 0
                extvbars[i1ep].SetForeRightBL(i1ep == 0, wvtrackpair.vbar)
                if wvtrackpair.i1 == 0:
                    assert extvbars[wvtrackpair.i1].cellmarkright is None
                    extvbars[wvtrackpair.i1].cellmarkright = self.izonecellmarkmap[wvtrackpair.izoneR]
                else:
                    assert extvbars[wvtrackpair.i1].cellmarkleft is None
                    extvbars[wvtrackpair.i1].cellmarkleft = self.izonecellmarkmap[wvtrackpair.izoneR]

        self.vbm.bars = [wvtrackpair.vbar  for wvtrackpair in wvtrackpairs] + extvbars
                
        # check we have all consistent cellmarks
        for bar in self.vbm.bars:
            assert DVBchecksetcellmark(bar, bar.nodeback)
            assert DVBchecksetcellmark(bar, bar.nodefore)



def DVBchecksetcellmark(bar, node):
    assert not bar.bbardeleted
    if bar.GetForeRightBL(bar.nodeback == node) is None:
        return True
    lbar, lnode = bar, node
    Dcellmark = lbar.GetCellMarkRightL(lnode == lbar.nodeback)
    Dcount = 0
    while True:
        lnode = lbar.GetNodeFore(lbar.nodeback == lnode)
        lbar = lbar.GetForeRightBL(lbar.nodefore == lnode)
        if Dcellmark != lbar.GetCellMarkRightL(lnode == lbar.nodeback):
            return False
        if lbar == bar:
            break
        Dcount += 1
        assert Dcount < 1000
    return True
            
        
def GetBarCycles(bm):
    extbarset = set(bar  for bar in bm.bars  if not bar.bbardeleted and bar.barbackleft is None or bar.barforeright is None)
    extbarcycles = [ ]
    while extbarset:
        bar = min(extbarset, key=lambda X:(X.nodeback.i, X.nodefore.i))  # avoid arbitraryness of pop()
        extbarset.remove(bar)
        node = bar.GetNodeFore(bar.barbackleft is not None)
        extbarcycle = [ ]
        sbar = bar
        while True:
            node = bar.GetNodeFore(bar.nodeback == node)
            extbarcycle.append((bar, node))
            while bar.GetForeRightBL(bar.nodefore == node) is not None:
                bar = bar.GetForeRightBL(bar.nodefore == node)
            if bar == sbar:
                break
            extbarset.remove(bar)
        extbarcycle.reverse()
        extbarcycles.append(extbarcycle)
    return extbarcycles


def MakeColpointNode(vbdrivegeometry, vbm, bar):
    assert not bar.bbardeleted
    if bar.cellmarkright is None or bar.cellmarkleft is None:
        return None
    izoneleft = bar.cellmarkleft.cellcolour
    izoneright = bar.cellmarkright.cellcolour
    
    isegleft = vbdrivegeometry.GetIseg(izoneleft)
    isegright = vbdrivegeometry.GetIseg(izoneright)

    if vbdrivegeometry.AreSegsAdjI(isegleft, isegright) != 0:
        return None
    
    vbsegleft = vbdrivegeometry.GetVBiseg(isegleft)
    vbsegright = vbdrivegeometry.GetVBiseg(isegright)
    if vbsegleft.IsLine() and vbsegright.IsLine():
        return None

    # bisection line between two corners
    if vbsegleft.IsZeroArc() and vbsegright.IsZeroArc():
        vp = vbsegright.pcen - vbsegleft.pcen
        r = vp.Len()*0.5
        mp = vbsegleft.pcen + vp*0.5
        vb = bar.nodefore.p - bar.nodeback.p
        lam = P3.Dot(P3.ConvertGZ(mp, 0.0) - bar.nodeback.p, vb) / vb.Lensq()
        if 0 < lam < 1:
            return NewNodeR(vbm, P3.ConvertGZ(mp, 0.0), r)
        return None
    
    # parabolic case
    vbline = vbsegleft if vbsegleft.IsLine() else vbsegright
    vbarc = vbsegright if vbsegleft.IsLine() else vbsegleft
    assert vbarc.IsZeroArc()

    avec = P2.CPerp(vbline.vperpnorm)
    parccen = P2(P2.Dot(avec, vbarc.pcen), P2.Dot(vbline.vperpnorm, vbarc.pcen) - vbline.vperpnormdot)
    pback = P2(P2.DotLZ(avec, bar.nodeback.p), P2.DotLZ(vbline.vperpnorm, bar.nodeback.p) - vbline.vperpnormdot)
    pfore = P2(P2.DotLZ(avec, bar.nodefore.p), P2.DotLZ(vbline.vperpnorm, bar.nodefore.p) - vbline.vperpnormdot)
    if (pback.u < parccen.u) != (pfore.u < parccen.u):
        r = parccen.v*0.5
        mp = avec * parccen.u + vbline.vperpnorm * (r + vbline.vperpnormdot)
        return NewNodeR(vbm, mp, abs(r))
    return None


def InsertColPointNodes(vbdrivegeometry, vbm):
    colnodes = [ ]
    for bar in vbm.bars:
        if bar.bbardeleted:
            continue
        cnode = MakeColpointNode(vbdrivegeometry, vbm, bar)
        if not cnode:
            continue
        colnodes.append((bar, cnode))
        
    for bar, cnode in colnodes:
         vbm.InsertNodeIntoBarF(bar, cnode, bcolinear=False)
