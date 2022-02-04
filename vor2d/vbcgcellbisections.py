from basicgeo import P2, Along
from .vbdrivegeometrycalculations import FindBisectVectorNorm, MakeCenpt3, MakeCenpt3Z

class CalcCutBisect2Obj:
    def __init__(self, barnode1, barnode2):
        bar1, node1 = barnode1
        bar2, node2 = barnode2
        self.izone = node2.izone
        assert self.izone == bar1.GetNodeFore(bar1.nodeback == node1).izone
        assert self.izone != node1.izone
        assert self.izone != bar2.GetNodeFore(bar2.nodeback == node2).izone

        self.bar1, self.node1 = bar1, node1
        self.bar2, self.node2 = bar2, node2
        assert self.bar1.nodemid is not None
        assert self.bar2.nodemid is not None
        # self.bisectvec1 = FindBisectVectorNorm
        # self.bisectvec2 = FindBisectVectorNorm
        
    def FindHalfPointInizone(self):
        return (self.bar1.nodemid.p + self.bar2.nodemid.p)*0.5
        
    # b3waycase thrown in to force ability to split a cell
    def CalcBisectorDivergence(self, vbs, vbbarcell, b3waycase):
        if self.bisectvec1 is None or self.bisectvec2 is None:
            return -1, 0, 0
            
        # do the two bisectvecs diverge or converge?
        pvm = self.bar2.nodemid.p - self.bar1.nodemid.p
        nmch = P2.Dot(pvm, P2.APerp(self.bisectvec1))
        nmchv = P2.Dot(self.bisectvec2, P2.APerp(self.bisectvec1))
        if (nmch > 0.0) != (nmchv > 0.0) and not b3waycase:
            return -1, 0, 0
            
        vc = P2.ZNorm(P2.APerp(self.bisectvec1 + self.bisectvec2))
        vch = P2.Dot(vc, (self.bar2.nodemid.p + self.bar1.nodemid.p)*0.5)
        
        vcdots = [ P2.Dot(vc, node.p) - vch  for bar, node in vbbarcell.barnodes ]
        divquality = min(-min(vcdots), max(vcdots))  # was abs(nmchv)

        return divquality, vc, vch
        
    def CalcBisectorJoin(self, vbs):
        # use vbs.FindBisectPoint on bar1 and bar2 to get the tangencies
        # and be able to tell if they are divergent or convergent to what point
        # to find if they should be split out-wise
        izoneback = self.node1.izone
        izonefore = self.bar2.GetNodeFore(self.bar2.nodeback == self.node2).izone
        self.cenpt, self.r = MakeCenpt3Z(vbs, (izoneback, self.izone, izonefore), True)  # speculating
        if self.cenpt and self.bisectvec1 and self.bisectvec2:
            self.pdist2 = P2.Dot(self.bisectvec2, self.cenpt - self.bar2.nodemid.p)
            self.pdist1 = P2.Dot(self.bisectvec1, self.cenpt - self.bar1.nodemid.p)
            self.bbackward = (self.pdist1 <= 0) or (self.pdist2 <= 0)
        else:
            self.bbackward = True
            
        
class CalcCutLine4Obj:
    def __init__(self, vbbarcell, vbs, b3waycase, Dnextvbarstack):
        self.Dnextvbarstack = Dnextvbarstack
        self.izonebars = [ CalcCutBisect2Obj(vbbarcell.barnodecuts[i], vbbarcell.barnodecuts[i+1  if i+1 != len(vbbarcell.barnodecuts)  else 0])  for i in range(len(vbbarcell.barnodecuts)) ]
        self.b3waycase = b3waycase
        
        # does one of the izones appear twice?
        izonecount = { }
        for bar, node in vbbarcell.barnodecuts:
            izonecount[node.izone] = izonecount.get(node.izone, 0) + 1
        izone2s = [ izone  for izone, c in izonecount.items()  if c >= 2 ]
        if izone2s:
            vc, vch = self.calcsplitthroughizone(izone2s[0], vbbarcell)
        else:
            vc, vch = self.calcsplitsomeizone(vbbarcell, vbs)
        self.ktnodetop, self.ktbartop, self.ktbartoplam, self.ktnodebot, self.ktbarbot, self.ktbarbotlam = vbbarcell.VBCellCutBar(vc, vch)

    # there's a (sometimes very narrow) vbpoly channel through the cell which we make a direct connection across
    def calcsplitthroughizone(self, izone2, vbbarcell):
        izonebars2 = [ izonebar  for izonebar in self.izonebars  if izonebar.izone == izone2 ]
        ptop = izonebars2[0].FindHalfPointInizone()
        pbot = izonebars2[1].FindHalfPointInizone()
        
        # joining line between the two points
        vc = P2.APerp(pbot - ptop)
        vch = P2.Dot(vc, ptop)
        return vc, vch

    def calcsplitsomeizone(self, vbbarcell, vbs):
        # fetch the pairs of bisectors on either side of each izone
        for i, izonebar in enumerate(self.izonebars):
            izonebar.bisectvec1 = FindBisectVectorNorm(vbs, izonebar.bar1, izonebar.node1)
            iprev = i-1  if i!=0  else len(self.izonebars)-1
            izonebarprev = self.izonebars[iprev]
            assert (izonebar.bar1, izonebar.node1) == (izonebarprev.bar2, izonebarprev.node2)
            izonebarprev.bisectvec2 = izonebar.bisectvec1
            
        # find all diverging cases neighbouring bisectors and find which is furthest from the sides
        izonebardivergent, izbd = max(((izonebar, izonebar.CalcBisectorDivergence(vbs, vbbarcell, self.b3waycase))  for izonebar in self.izonebars), key=lambda X:X[1][0])

        if izbd[0] > 0.0 or self.b3waycase:
            Dvc, Dvch = P2.ZNorm(izbd[1]), izbd[2]/izbd[1].Len()
            return izbd[1], izbd[2]  # vc, vch

        assert len(vbbarcell.barnodecuts) == len(self.izonebars)
        for izonebar in self.izonebars:
            izonebar.CalcBisectorJoin(vbs)
            
        # look for most distant opposing zones starting with the known merging case
        i0 = min([i  for i in range(len(self.izonebars))  if not self.izonebars[i].bbackward], default=None, key=lambda X:self.izonebars[X].r)
        if i0 is not None:
            i0neigh = [i0-1 if i0 != 0 else len(self.izonebars)-1, i0, i0+1 if i0 != len(self.izonebars)-1 else 0]
            i1 = min([i  for i in range(len(self.izonebars))  if i not in i0neigh and not self.izonebars[i].bbackward], default=None, key=lambda X:self.izonebars[X].r)
            if i1 is not None:
                cenpt0 = self.izonebars[i0].cenpt
                cenpt1 = self.izonebars[i1].cenpt
                vc = P2.ZNorm(cenpt1 - cenpt0)
                vch = P2.Dot(vc, Along(0.5, cenpt0, cenpt1))
                return vc, vch

        # find two subsequent Cenpt3 positions and bisect them 
        assert len(vbbarcell.barnodecuts) >= 4
        def MakeHypCenpt3(i, iknockout):
            barnodecuts3 = [ ]
            im1 = i-1 if i != 0 else len(vbbarcell.barnodecuts)-1
            ip1 = i+1 if i != len(vbbarcell.barnodecuts)-1 else 0
            if iknockout == -1:
                im2 = im1-1 if im1 != 0 else len(vbbarcell.barnodecuts)-1
                barnodecuts3.append(vbbarcell.barnodecuts[im2])
            barnodecuts3.append(vbbarcell.barnodecuts[im1])
            if iknockout == 0:
                barnodecuts3.append(vbbarcell.barnodecuts[i])
            barnodecuts3.append(vbbarcell.barnodecuts[ip1])
            if iknockout == 1:
                ip2 = ip1+1 if ip1 != len(vbbarcell.barnodecuts)-1 else 0
                barnodecuts3.append(vbbarcell.barnodecuts[ip2])
            assert len(barnodecuts3) == 3
            cenpt, r = MakeCenpt3(vbs, barnodecuts3, True)   # speculating
            return ((i, iknockout), cenpt, r)
            
        minr = min((bar.nodemid.r  for bar, node in vbbarcell.barnodecuts))
        minr = 0.0  # override
        
        pp3cenpts = [ MakeHypCenpt3(i, 0)  for i in range(len(vbbarcell.barnodecuts)) ]
        (i, knockout), cenpt, r = min((pp3  for pp3 in pp3cenpts  if pp3[2] >= minr), key=lambda X:X[2])
        pp3cenpts.append(MakeHypCenpt3(i, -1))
        pp3cenpts.append(MakeHypCenpt3(i, +1))
        i1, cenpt1, r1 = min((pp3  for pp3 in pp3cenpts  if pp3[2] >= minr and pp3[0] != (i, 0)), key=lambda X:X[2])
        if cenpt1 is None or cenpt is None:
            print("bark5", pp3cenpts, minr)
            print(min(pp3  for pp3 in pp3cenpts  if pp3[2] >= minr))
            sendactivity(contours=[[node.p  for bar, node in vbbarcell.barnodes]])
            sendactivity(points=[cenpt  for ((i, iknockout), cenpt, r) in pp3cenpts  if cenpt is not None])
            sendactivity(contours=[[bar.nodemid.p  for bar, node in vbbarcell.barnodecuts]], materialnumber=3)
            print([bar.nodemid.r  for bar, node in vbbarcell.barnodecuts])
            self.bari = None
            assert False
            return

        # bisection line between the two points
        vc = P2.ZNorm(cenpt1 - cenpt)
        vch = P2.Dot(vc, Along(0.5, cenpt, cenpt1))
        return vc, vch

        
    def makesplitnodes(self, btop, bm, lamendgap):
        if btop:
            ktnode, ktbar, ktbarlam = self.ktnodetop, self.ktbartop, self.ktbartoplam
        else:
            ktnode, ktbar, ktbarlam = self.ktnodebot, self.ktbarbot, self.ktbarbotlam
        assert not ktbar.bbardeleted
        ktnodefore = ktbar.GetNodeFore(ktbar.nodeback == ktnode)
        
        if lamendgap < ktbarlam < 1 - lamendgap:
            splitnode = bm.InsertNodeIntoBarF(ktbar, bm.NewNode(Along(ktbarlam, ktnode.p, ktnodefore.p)), True)
            assert ktbar.bbardeleted
            leadsplitbar = ktbar.GetForeRightBL(ktnode == ktbar.nodefore)
            assert not leadsplitbar.bbardeleted
            bnewnode = True
            
        else:
            #print("applyinglamendgap", ktbarlam, lamendgap, Along(ktbarlam, ktnode.p, ktnodefore.p))
            if ktbarlam < 0.5:
                splitnode = ktnode
                if ktnode == ktbar.nodeback:
                    leadsplitbar = ktbar.GetBarBackRight()
                else:
                    leadsplitbar = ktbar.GetBarForeLeft()
                assert splitnode == leadsplitbar.nodeback or splitnode == leadsplitbar.nodefore
            else:
                splitnode = ktnodefore
                leadsplitbar = ktbar
                assert splitnode == leadsplitbar.nodeback or splitnode == leadsplitbar.nodefore
            bnewnode = False

        if btop:
            self.leadsplitbartop, self.splitnodetop = leadsplitbar, splitnode
            assert self.splitnodetop == self.leadsplitbartop.nodeback or self.splitnodetop == self.leadsplitbartop.nodefore
            self.bsplitnnodetopnew = bnewnode
        else:
            self.leadsplitbarbot, self.splitnodebot = leadsplitbar, splitnode            
            assert self.splitnodebot == self.leadsplitbarbot.nodeback or self.splitnodebot == self.leadsplitbarbot.nodefore
            self.bsplitnnodebotnew = bnewnode
        return bnewnode
            
    def makecellsplittingbar(self, bm):
        if bm.DTestColinearityF(self.splitnodetop, self.leadsplitbartop, self.splitnodebot, self.leadsplitbarbot) or bm.DTestColinearityF(self.splitnodebot, self.leadsplitbarbot, self.splitnodetop, self.leadsplitbartop):
            return None
        if self.splitnodetop.i < self.splitnodebot.i:
            splitbar = bm.MakeBarBetweenNodesF(self.splitnodetop, self.leadsplitbartop, self.splitnodebot, self.leadsplitbarbot)
        else:
            splitbar = bm.MakeBarBetweenNodesF(self.splitnodebot, self.leadsplitbarbot, self.splitnodetop, self.leadsplitbartop)
        return splitbar
        
    def makesplitbarsandnodes(self, bm, lamendgap):
        newsplitnodes, newsplitbars = [ ], [ ]
        self.makesplitnodes(True, bm, lamendgap)
        self.makesplitnodes(False, bm, lamendgap)
    
        if self.bsplitnnodetopnew:
            newsplitnodes.append(self.splitnodetop)
            newsplitbars.append(self.leadsplitbartop)
            newsplitbars.append(self.leadsplitbartop.GetForeRightBL(self.leadsplitbartop.nodefore == self.splitnodetop))
        else:
            newsplitbars.append(self.ktbartop)
            
        if self.bsplitnnodebotnew:
            newsplitnodes.append(self.splitnodebot)
            newsplitbars.append(self.leadsplitbarbot)
            newsplitbars.append(self.leadsplitbarbot.GetForeRightBL(self.leadsplitbarbot.nodefore == self.splitnodebot))
        else:
            newsplitbars.append(self.ktbarbot)
        return newsplitnodes, newsplitbars
