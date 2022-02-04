from basicgeo import P2, Along, I1
from .vbconstructionbarmesh import SetCellmark, CellMark, Node
from .vbdrivegeometrycalculations import FindBisectVectorNorm, FindBisectPoint, MakeCenpt3, MakeCenpt3Z, DetectMissingThreeway
from .vbcgcellbisections import CalcCutBisect2Obj, CalcCutLine4Obj
import time

class VBbarcell:
    def __init__(self, bar, node):
        assert not bar.bbardeleted
        assert bar.GetForeRightBL(bar.nodeback == node) is not None  # doesn't go off edge
        self.bar = bar
        self.node = node
        #self.barnodes = [ ]
        #self.barnodecuts = [ ]
        #self.Durg, self.Dvrg
        
    def GetCellMark(self):
        return self.bar.GetCellMarkRightL(self.node == self.bar.nodeback)

    def MakeBarNodeCutsOfCell(self):
        assert not self.bar.bbardeleted
        self.barnodes = [ ]
        self.barnodecuts = [ ]
        Dcellmark = self.GetCellMark()
        bar, node = self.bar, self.node
        Dcount = 0
        while True:
            self.barnodes.append((bar, node))
            if (bar.nodeback.izone != bar.nodefore.izone):
                if self.barnodecuts:
                    assert node.izone != self.barnodecuts[-1][1].izone
                    assert node.izone == self.barnodecuts[-1][0].GetNodeFore(self.barnodecuts[-1][0].nodeback == self.barnodecuts[-1][1]).izone
                self.barnodecuts.append((bar, node))
            node = bar.GetNodeFore(bar.nodeback == node)
            bar = bar.GetForeRightBL(bar.nodefore == node)
            assert bar is not None and not bar.bbardeleted
            if bar == self.bar:
                assert node == self.node
                break
            assert bar.GetCellMarkRightL(node == bar.nodeback) == Dcellmark
            Dcount += 1
            assert Dcount < 1000
            
        self.Durg = I1(min(node.p[0]  for bar, node in self.barnodes), max(node.p[0]  for bar, node in self.barnodes))
        self.Dvrg = I1(min(node.p[1]  for bar, node in self.barnodes), max(node.p[1]  for bar, node in self.barnodes))
    
    def CheckIsolationSegment(self, vbsegmentslist):
        #for 
        #        def GetIseg(self, izone):
        #return (izone % len(self.vbsegmentslist), izone // len(self.vbsegmentslist))
        pass
        
    # returns equivalent of two calls to PolyCutBar
    def VBCellCutBar(self, vc, vch):
        vcdots = [ P2.Dot(vc, node.p) - vch  for bar, node in self.barnodes ]
        if min(vcdots) > -1e-4 or max(vcdots) < 1e-4:
            print("bark4", min(vcdots), max(vcdots), vcdots)
            #sendactivity(contours=[[node.p  for bar, node in self.barnodes]], materialnumber=1)
            Dvc, Dvch = P2.ZNorm(vc), vch*vc.Len()
            #sendactivity(contours=[[Dvc*Dvch + P2.APerp(Dvc)*100, Dvc*(Dvch/Dvc.Lensq()) - P2.APerp(Dvc)*100]], materialnumber=1)
            #ddd
            dvch = (max(vcdots) + min(vcdots))*0.5
            vch += dvch
            vcdots = [ P2.Dot(vc, node.p) - vch  for bar, node in self.barnodes ]
        
        for i in range(len(self.barnodes)):
            ip = i+1  if i+1 != len(self.barnodes)  else 0
            if vcdots[i] < 0 <= vcdots[ip]:
                ktnodetop, ktbartop, ktbartoplam = self.barnodes[i][1], self.barnodes[i][0], -vcdots[i] / (vcdots[ip] - vcdots[i])
            if vcdots[i] >= 0 > vcdots[ip]:
                ktnodebot, ktbottop, ktbarbotlam = self.barnodes[i][1], self.barnodes[i][0], -vcdots[i] / (vcdots[ip] - vcdots[i])
        return ktnodetop, ktbartop, ktbartoplam, ktnodebot, ktbottop, ktbarbotlam
        
    def DInteriorness(self, p):
        res = -1
        winding = 0
        nodeprev = self.barnodes[-1][1]
        for bar, node in self.barnodes:
            ndist = (node.p - p).Len()
            if res == -1 or ndist < res:
                res = ndist
            v = node.p - nodeprev.p
            lam = P2.Dot(p - nodeprev.p, v) / v.Lensq()
            if 0 < lam < 1:
                bdist = (Along(lam, nodeprev.p, node.p) - p).Len()
                if bdist < res:
                    res = bdist
            if (node.p.v < p.v) != (nodeprev.p.v < p.v):
                vlam = (p.v - nodeprev.p.v)/(node.p.v - nodeprev.p.v)
                npu = Along(vlam, nodeprev.p.u, node.p.u)
                if npu > p.u:
                    winding += 1 if nodeprev.p.v>node.p.v else -1
            nodeprev = node
        assert winding == 0 or winding == 1
        assert res >= 0
        return res * (-1 if winding == 0 else 1)


class VBconstructiongeometry:
    def __init__(self, lvbdrivegeometry, lbm):
        self.vbdrivegeometry = lvbdrivegeometry
        self.bm = lbm
        self.Dnextvbarstack = 0
        self.vbbarcells = [ ]

        self.vbsclosesttime = 0
        self.vbstotaltime = 0
        self.vbsclosestcalls = 0
        self.lamendgap = 0.0001
        self.Dtooshortbars = [ ]
        self.shortestsubdivbarlength = 1e-9
        
    def goconstruct(self, maxDnextvbarstack=-1):
        totalctime = time.clock()
        if self.Dnextvbarstack == 0:
            ctime = time.clock()
            self.vbdrivegeometry.MakePointZoneVBS(self.bm.nodes)
            self.vbsclosesttime += time.clock() - ctime
            self.vbsclosestcalls += len(self.bm.nodes)
            
            self.barstack = [ bar  for bar in self.bm.bars  if not bar.bbardeleted and bar.nodeback.izone != bar.nodefore.izone ]
            self.bm.maxcellcolour = -1

        while self.barstack and (maxDnextvbarstack == -1 or self.Dnextvbarstack < maxDnextvbarstack):
            self.nextvbarstack()
        self.vbstotaltime += time.clock() - totalctime

    def InvalidateCellMarksOnBar(self, bar):
        oldcellmarkRight = bar.GetCellMarkRightL(True)
        oldcellmarkLeft = bar.GetCellMarkRightL(False)
        if bar.GetForeRightBL(True) is not None and oldcellmarkRight is not None:
            assert self.vbbarcells[oldcellmarkRight.cellcolour] is not None
            self.vbbarcells[oldcellmarkRight.cellcolour] = None
            SetCellmark(bar, bar.nodeback, None, oldcellmarkRight)
        if bar.GetForeRightBL(False) is not None and bar.GetCellMarkRightL(False) is not None:
            assert self.vbbarcells[oldcellmarkLeft.cellcolour] is not None
            self.vbbarcells[oldcellmarkLeft.cellcolour] = None
            SetCellmark(bar, bar.nodefore, None, oldcellmarkLeft)

    def MakeNodeMid(self, bar):
        pb, biscontained = FindBisectPoint(self.vbdrivegeometry, bar.nodeback, bar.nodefore)
        if not pb:
            pb = (bar.nodeback.p + bar.nodefore.p)*0.5
        elif not biscontained:
            v = bar.nodefore.p + bar.nodeback.p
            lam = P2.Dot(pb - bar.nodeback.p, v)/v.Lensq()
            if not (0.1 < lam < 0.9):
                pb = (bar.nodeback.p + bar.nodefore.p)*0.5
        else:
            assert pb != bar.nodeback.p and pb != bar.nodefore.p, (biscontained, bar.nodemid and bar.nodemid.p)
            
        assert pb != bar.nodeback.p and pb != bar.nodefore.p, (pb, bar.nodeback.p, bar.nodefore.p)
        ctime = time.clock()
        pzm = Node(pb, -1)
        self.vbdrivegeometry.VBclosest(pzm, pb)
        self.vbsclosesttime += time.clock() - ctime
        self.vbsclosestcalls += 1
        
        # force a better subdivision point that is not in either side
        if not biscontained and (pzm.izone == bar.nodeback.izone or pzm.izone == bar.nodefore.izone):
            p0, izone0 = bar.nodeback.p, bar.nodeback.izone
            p1, izone1 = bar.nodefore.p, bar.nodefore.izone
            Nnewpzms = 0
            while (p0 - p1).Len() > 1e-9:
                if pzm.izone == izone0:
                    p0 = pb
                elif pzm.izone == izone1:
                    p1 = pb
                else:
                    break
                pb = (p0 + p1)*0.5
                ctime = time.clock()
                pzm = Node(pb, -1)
                self.vbdrivegeometry.VBclosest(pzm, pb)
                self.vbsclosesttime += time.clock() - ctime
                self.vbsclosestcalls += 1
                Nnewpzms += 1
            if Nnewpzms>5:
                print("Nnewpzms", Nnewpzms, self.vbdrivegeometry.GetIseg(izone0), self.vbdrivegeometry.GetIseg(pzm.izone), self.vbdrivegeometry.GetIseg(izone1), pzm.r)
            
        bar.nodemid = self.bm.NewNode(pb)
        bar.nodemid.izone, bar.nodemid.r = pzm.izone, pzm.r
        return biscontained
        
    def nextvbarstack(self):
        bar = self.barstack.pop()
        self.Dnextvbarstack += 1
        
        if bar.bbardeleted:  
            return
            
        # if non-splitting ones are in the list it's for the reason that we need to check the vbbarcells on either side again anyway
        # if bar.nodeback.izone == bar.nodefore.izone:  
        #    return
        
        if bar.nodeback.izone != bar.nodefore.izone and not bar.nodemid:
            biscontained = self.MakeNodeMid(bar)

            # nodemid is not the final boundary between the two endpoints of the bar
            if not biscontained or (bar.nodemid.izone != bar.nodeback.izone and bar.nodemid.izone != bar.nodefore.izone):
                assert bar.nodemid.p != bar.nodeback.p and bar.nodemid.p != bar.nodefore.p, bar.nodemid.p
                if (bar.nodeback.p - bar.nodefore.p).Len() < self.shortestsubdivbarlength:
                    print("skipping too short; what to do?", bar.nodeback.p)
                    self.Dbarshort = bar
                    self.Dtooshortbars.append(bar)
                    #assert False
                    return

                self.InvalidateCellMarksOnBar(bar)

                # drop the new split bars back into the stack for further processing
                self.Dbarnodemid = bar
                self.bm.InsertNodeIntoBarF(bar, bar.nodemid, True)
                assert bar.bbardeleted
                self.barstack.append(bar.barforeright)
                self.barstack.append(bar.barbackleft)
                return
            else:
                pass  # drop through into having a no further subdividing bar we should inspect the cells of
            
        assert not bar.bbardeleted
        assert ((bar.nodeback.izone == bar.nodefore.izone) or bar.nodemid is not None)

        # create the barcells on either side of the bar (will bail out of anything happens)
        vbbarcells2 = [ ]
        if bar.GetForeRightBL(True) is not None:
            vbbarcells2.append(VBbarcell(bar, bar.nodeback))
        if bar.GetForeRightBL(False) is not None:
            vbbarcells2.append(VBbarcell(bar, bar.nodefore))
        
        fcl4 = None  # drop through if any subdividing happens on one or other side
        while vbbarcells2:
            vbbarcell = vbbarcells2.pop()
            if vbbarcell.GetCellMark() is not None:
                continue    # already done
            vbbarcell.MakeBarNodeCutsOfCell()  
            self.Dvbbarcell = vbbarcell
            if len(vbbarcell.barnodecuts) == 0:
                continue    # no cuts here
            if sum(1  for bar, node in vbbarcell.barnodecuts  if bar.nodemid is None) != 0:
                continue    # not all sides have been processed to the limit of subdivision.  once we rotate through to the final one we can do it
                
            assert len(vbbarcell.barnodecuts) >= 2
            if len(vbbarcell.barnodecuts) >= 4:
                fcl4 = CalcCutLine4Obj(vbbarcell, self.vbdrivegeometry, False, self.Dnextvbarstack)
                break  # drop out to work on subdividing this
                
            assert len(vbbarcell.barnodecuts) >= 2
            if len(vbbarcell.barnodecuts) == 3: 
                vbbarcell.cenpt, vbbarcell.cenptR = MakeCenpt3(self.vbdrivegeometry, vbbarcell.barnodecuts, False)   # real case
                if vbbarcell.cenptR < 0:  
                    if DetectMissingThreeway(self.vbdrivegeometry, vbbarcell):
                        print("DetectedMissingThreeway")
                    else:
                        print("proj failure", vbbarcell.cenpt, vbbarcell.cenptR)
                        assert False
                    fcl4 = CalcCutLine4Obj(vbbarcell, self.vbdrivegeometry, True, self.Dnextvbarstack)
                    break  # drop out to work on subdividing this
                assert vbbarcell.Durg.Contains(vbbarcell.cenpt.u) and vbbarcell.Dvrg.Contains(vbbarcell.cenpt.v), (vbbarcell.cenpt, vbbarcell.Durg, vbbarcell.Dvrg)
                assert vbbarcell.DInteriorness(vbbarcell.cenpt) >= -1e-5, vbbarcell.cenpt
                
            # put new cell mark onto the face and corresponding vbbarcell into the indexed list
            self.bm.maxcellcolour += 1
            newcellmark = CellMark(self.bm.maxcellcolour)
            SetCellmark(vbbarcell.bar, vbbarcell.node, newcellmark, None)
            assert len(self.vbbarcells) == newcellmark.cellcolour
            self.vbbarcells.append(vbbarcell)
            assert self.vbbarcells[newcellmark.cellcolour].GetCellMark() == newcellmark
            
        if fcl4 is None:
            return
            
        # cut and insert a bar into this cell
        self.InvalidateCellMarksOnBar(fcl4.ktbartop)
        self.InvalidateCellMarksOnBar(fcl4.ktbarbot)

        newsplitnodes, newsplitbars = fcl4.makesplitbarsandnodes(self.bm, lamendgap=self.lamendgap)

        ctime = time.clock()
        self.vbdrivegeometry.MakePointZoneVBS(newsplitnodes)
        self.vbsclosesttime += time.clock() - ctime
        self.vbsclosestcalls += len(newsplitnodes)

        splitbar = fcl4.makecellsplittingbar(self.bm)
        if splitbar is None:
            print("ffail!!!!!!")   # should be an assert and make us try it again, maybe in some better directions
            return  # fail
            
        # push all the new stuff back onto the stack
        self.barstack.append(splitbar)
        if not vbbarcell.bar.bbardeleted:
            self.barstack.append(vbbarcell.bar)
        self.barstack.extend(newsplitbars)

    def PlotBarMesh(self, sendactivity, bpoints=True):
        def DGetMidP(bar):
            return bar.nodemid.p if bar.nodemid else (bar.nodefore.p + bar.nodeback.p)*0.5
        if bpoints:
            vds = { }
            for nd in self.bm.nodes:
                izone = nd.izone
                if izone not in vds:
                    vds[izone] = [ ]
                vds[izone].append(P2(nd.p.x, nd.p.y))
            for i, l in vds.items():
                sendactivity("points", points=l, materialnumber=(i % 5))
            sendactivity("points", points=[bar.nodemid.p  for bar in self.bm.bars  if (not bar.bbardeleted) and bar.nodemid ])
        
        conts = [ ]
        mconts = [ ]
        for vbbarcell in self.vbbarcells:
            if vbbarcell is None:
                continue
            if len(vbbarcell.barnodecuts) == 2:
                conts.append((DGetMidP(vbbarcell.barnodecuts[0][0]), DGetMidP(vbbarcell.barnodecuts[1][0])))
            elif len(vbbarcell.barnodecuts) == 3 and vbbarcell.cenpt is not None:
                for bar, node in vbbarcell.barnodecuts:
                    #if bar.nodemid:
                        conts.append((bar.nodemid.p, vbbarcell.cenpt))
            else:
                pbar, pnode = vbbarcell.barnodecuts[-1]
                avgpt = sum((DGetMidP(bar)  for bar, node in vbbarcell.barnodecuts), P2(0,0)) * (1/len(vbbarcell.barnodecuts))
                for bar, node in vbbarcell.barnodecuts:
                    if bar.nodemid:
                        mconts.append((bar.nodemid.p, avgpt))
                
                #for i, (bar, node) in enumerate(vbbarcell.barnodecuts):
                #    if not pbar.nodemid or not bar.nodemid:  continue  # hack to handle bad insertions done temporarily
                #    conts.append((pbar.nodemid.p, bar.nodemid.p))
                #    pbar, pnode = bar, node
        sendactivity(contours=conts, materialnumber=3)
        sendactivity(contours=mconts, materialnumber=1)

