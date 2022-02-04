from barmesh import PointZone

# PointZone { izone: {PZ_WITHIN_R, PZ_BEYOND_R}, r, v }

class ImpingementShape:
    def __init__(self, tbarmesh, tboxing, DPZshape, DLPZshape):
        self.tbarmesh = tbarmesh
        self.tboxing = tboxing
        self.DPZshape = DPZshape
        self.DLPZshape = DLPZshape
        self.hitreg = [0]*len(tbarmesh.bars)
        self.nhitreg = 0

    def Isb2dcontournormals(self):
        return False
    
    def DistP(self, pz, p):
        dpz = self.DPZshape(p, pz.r, self)
        
        for ixy in self.tboxing.CloseBoxeGenerator(p.x, p.x, p.y, p.y, pz.r):
            tbox = self.tboxing.GetTriangleBox(ixy)
            #for i in self.tboxing.SlicePointisZ(tbox.pointis, p.z-pz.r, p.z+pz.r):
            for i in tbox.pointis:
                dpz.DistPpointPZ(self.tboxing.GetNodePoint(i))
                
        self.nhitreg += 1
        for ixy in self.tboxing.CloseBoxeGenerator(p.x, p.x, p.y, p.y, pz.r):
            tbox = self.tboxing.GetTriangleBox(ixy)
            for i in tbox.edgeis:
                if self.hitreg[i] != self.nhitreg:
                    dpz.DistPedgePZ(*self.tboxing.GetBarPoints(i))
                    self.hitreg[i] = self.nhitreg
                    
        # then do the triangles
        self.nhitreg += 1
        for ixy in self.tboxing.CloseBoxeGenerator(p.x, p.x, p.y, p.y, pz.r):
            tbox = self.tboxing.GetTriangleBox(ixy)
            for i in tbox.triangleis:
                if self.hitreg[i] != self.nhitreg:
                    dpz.DistPtrianglePZ(*self.tboxing.GetTriPoints(i))
                    self.hitreg[i] = self.nhitreg
                    
        pz.r = dpz.r
        pz.v = dpz.v
        

    def Cutpos(self, p, vp, cp, r):  # point, vector, cp=known close point to narrow down the cutoff search
        dlpz = self.DLPZshape(p, vp, r, self)  
        if cp is not None:
            dlpz.DistLamPpointPZ(cp)
            assert dlpz.lam != 2.0
        
        # solve |p + vp * lam - p| == r where Dist(p0) >= r >= Dist(p0 + vp)
        rexp = r + 0.01
        for ixy in self.tboxing.CloseBoxeGenerator(min(p.x, p.x+vp.x), max(p.x, p.x+vp.x), min(p.y, p.y+vp.y), max(p.y, p.y+vp.y), rexp):
            tbox = self.tboxing.GetTriangleBox(ixy)
            for i in self.tboxing.SlicePointisZ(tbox.pointis, min(p.z,p.z+vp.z)-rexp, max(p.z,p.z+vp.z)+rexp):
            #for i in tbox.pointis:
                dlpz.DistLamPpointPZ(self.tboxing.GetNodePoint(i))
                
        self.nhitreg += 1
        for ixy in self.tboxing.CloseBoxeGenerator(min(p.x, p.x+vp.x*dlpz.lam), max(p.x, p.x+vp.x*dlpz.lam), min(p.y, p.y+vp.y*dlpz.lam), max(p.y, p.y+vp.y*dlpz.lam), r + 0.01):
            tbox = self.tboxing.GetTriangleBox(ixy)
            for i in tbox.edgeis:
                if self.hitreg[i] != self.nhitreg:
                    dlpz.DistLamPedgePZ(*self.tboxing.GetBarPoints(i))
                    self.hitreg[i] = self.nhitreg
        
        self.nhitreg += 1
        for ixy in self.tboxing.CloseBoxeGenerator(min(p.x, p.x+vp.x*dlpz.lam), max(p.x, p.x+vp.x*dlpz.lam), min(p.y, p.y+vp.y*dlpz.lam), max(p.y, p.y+vp.y*dlpz.lam), r + 0.01):
            tbox = self.tboxing.GetTriangleBox(ixy)
            for i in tbox.triangleis:
                if self.hitreg[i] != self.nhitreg:
                    dlpz.DistLamPtrianglePZ(*self.tboxing.GetTriPoints(i))
                    self.hitreg[i] = self.nhitreg
        
        if __debug__:
            if not (dlpz.lam == 0.0 or dlpz.lam == 2.0):
                pz = PointZone(0, r + 1.1, None)
                self.DistP(pz, dlpz.p + dlpz.vp*dlpz.lam)
                assert abs(pz.r - r) < 0.002, ("cutposbad", pz.r, dlpz.lam, dlpz.p, dlpz.vp)
            
        return dlpz.lam


