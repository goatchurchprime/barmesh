from basicgeo import P3, P2
import math

def AlongP3Z(z, p0, p1):
    assert (p0.z < z) != (p1.z < z) or p0.z == z, (p0.z, z, p1.z)
    lam = (z - p0.z) / (p1.z - p0.z)
    assert 0.0 <= lam <= 1.0
    return P3(p0.x * (1 - lam) + p1.x * lam, p0.y * (1 - lam) + p1.y * lam, z)

def PedgeSlice(p0, p1, zlo, zhi):
    assert p0.z <= p1.z
    if p0.z < zlo:
        assert p1.z >= zlo
        p0z = AlongP3Z(zlo, p0, p1)
        p0zsliced = True
        if p1.z <= zhi:
            p1z = p1
            p1zsliced = False
        else:
            p1z = AlongP3Z(zhi, p0, p1)
            p1zsliced = True
    else:
        assert p0.z <= zhi
        assert p1.z >= zlo
        p0z = p0
        p0zsliced = False
        if p1.z <= zhi:
            p1z = p1
            p1zsliced = False
        else:
            p1z = AlongP3Z(zhi, p0, p1)
            p1zsliced = True
    assert zlo <= p0z.z <= zhi
    assert zlo <= p1z.z <= zhi
    return p0z, p0zsliced, p1z, p1zsliced

def PtriSlice(p0, p1, p2, z):
    assert p0.z <= p1.z <= p2.z
    assert p0.z <= z <= p2.z
    assert p0.z != p2.z, (p0.z, p1.z, p2.z, z)
    p0z = AlongP3Z(z, p0, p2)
    if p1.z <= z:
        p1z = AlongP3Z(z, p1, p2)
    else:
        p1z = AlongP3Z(z, p0, p1)
    return p0z, p1z
    

class DistPZC:
    def __init__(self, p, r, zlo, zhi):
        self.p = p
        self.r = r
        self.zlo = zlo
        self.zhi = zhi
        self.v = None
        self.Dpp = None
        
    def DistPpointPZ(self, p0):
        if not (self.zlo <= p0.z <= self.zhi):
            return
        lv = p0 - self.p
        lvf = P2(lv.x, lv.y)
        lvflen = lvf.Len()
        if lvflen < self.r:
            self.r = lvflen
            self.v = lv
            assert self.zlo <= self.p.z + self.v.z <= self.zhi, (self.zlo, self.p.z + self.v.z, self.zhi)
            self.Dpp = [p0]

    def DistPedgePZF(self, p0z, p1z):
        assert self.zlo <= p0z.z <= p1z.z <= self.zhi
        v = p1z - p0z
        vf = P2(v.x, v.y)
        vfsq = vf.Lensq()
        if vfsq == 0.0:
            return
        vflen = vf.Len()
        lv = self.p - p0z 
        dp = abs(P2.DotLZ(P2.CPerp(vf), lv) / vflen)
        if dp >= self.r:
            return
        lam = P2.DotLZ(vf, lv) / vfsq
        if 0 < lam < 1:
            self.r = dp
            self.v = v * lam - lv
            assert abs(P2.DotLZ(vf, self.v)) < 0.001
            assert abs(P2(self.v.x, self.v.y).Len() - self.r) < 0.001
            assert self.zlo <= self.p.z + self.v.z <= self.zhi, (self.zlo, self.p.z, self.p.z + self.v.z, self.zhi)
            self.Dpp = [p0z, p1z]

    def DistPedgePZ(self, p0, p1):
        if p0.z < self.zlo and p1.z < self.zlo:
            return
        if p0.z > self.zhi and p1.z > self.zhi:
            return
        if p0.z <= p1.z:
            p0z, p0zsliced, p1z, p1zsliced = PedgeSlice(p0, p1, self.zlo, self.zhi)
        else:
            p0z, p0zsliced, p1z, p1zsliced = PedgeSlice(p1, p0, self.zlo, self.zhi)
            
        if p0zsliced:
            self.DistPpointPZ(p0z)
        if p1zsliced:
            self.DistPpointPZ(p1z)
        self.DistPedgePZF(p0z, p1z)
                
    def DistPtrianglePZO(self, p0, p1, p2):
        assert p0.z <= p1.z <= p2.z
        if p0.z <= self.zlo <= p2.z:
            self.DistPedgePZF(*PtriSlice(p0, p1, p2, self.zlo))
        if p0.z <= self.zhi <= p2.z:
            self.DistPedgePZF(*PtriSlice(p0, p1, p2, self.zhi))

    def DistPtrianglePZ(self, p0, p1, p2):
        if p0.z <= p1.z:
            if p0.z <= p2.z:
                if p1.z <= p2.z:
                    self.DistPtrianglePZO(p0, p1, p2)
                else:
                    self.DistPtrianglePZO(p0, p2, p1)
            else:
                self.DistPtrianglePZO(p2, p0, p1)
        elif p0.z <= p2.z:
            self.DistPtrianglePZO(p1, p0, p2)
        elif p1.z <= p2.z:
            self.DistPtrianglePZO(p1, p2, p0)
        else:
            self.DistPtrianglePZO(p2, p1, p0)

class DistLamPZC:
    def __init__(self, p, vp, r, zlo, zhi):
        self.p = p
        self.vp = vp
        self.vpf = P2(vp.x, vp.y)
        self.vpfsq = self.vpf.Lensq()
        self.vpflen = self.vpf.Len()
        self.vpfperpnorm = P2.APerp(self.vpf) * (1.0 / self.vpflen)
        self.r = r
        self.rsq = r*r
        self.zlo = zlo
        self.zhi = zhi
        self.lam = 2.0
        self.Dllist = [ ]
        
    def DistLamPpointPZ(self, p0):
        assert self.zlo <= p0.z <= self.zhi, (self.zlo, p0.z, self.zhi)
        lv = self.p - p0
        sd = P2.DotLZ(self.vpfperpnorm, lv)
        sasq = self.rsq - sd*sd
        
        lvf = P2(p0.x - self.p.x, p0.y - self.p.y)
        if sasq < 0.0:
            return
        sa = math.sqrt(sasq / self.vpfsq)
        lamc = P2.Dot(self.vpf, lvf) / self.vpfsq
        if lamc + sa < 0.0:
            return
        laml = lamc - sa
        if laml < 0.0:
            self.lam = 0.0  # shouldn't happen
        elif laml < self.lam:
            self.lam = laml
            
    def DistLamPedgePZF(self, p0z, p1z):
        self.Dllist.append([p0z, p1z])
        assert self.zlo <= p0z.z <= p1z.z <= self.zhi
        v = p1z - p0z 
        vf = P2(v.x, v.y)
        vfsq = vf.Lensq()
        if vfsq == 0.0:
            return
        vflen = vf.Len()
        vfperpnorm = P2.CPerp(vf) * (1.0 / vflen)
        
        lv = self.p - p0z 
        
        dp = P2.DotLZ(vfperpnorm, lv)
        dpm = P2.DotLZ(vfperpnorm, self.vp)
        if dpm > 0.0:
            laml = (-self.r - dp) / dpm
        elif dpm < 0.0:
            laml = (self.r - dp) / dpm
        else:
            return
        if 0 < laml < self.lam:
            lvp = lv + self.vp * laml
            assert abs(abs(P2.DotLZ(vfperpnorm, lvp)) - self.r) < 0.001
            mu = P2.DotLZ(vf, lvp) / vfsq
            if 0 < mu < 1:
                assert abs(P2.DotLZ(vf, lvp - v * mu)) < 0.001
                self.lam = laml


    def DistLamPedgePZ(self, p0, p1):
        if p0.z < self.zlo and p1.z < self.zlo:
            return
        if p0.z > self.zhi and p1.z > self.zhi:
            return
        if p0.z <= p1.z:
            p0z, p0zsliced, p1z, p1zsliced = PedgeSlice(p0, p1, self.zlo, self.zhi)
        else:
            p0z, p0zsliced, p1z, p1zsliced = PedgeSlice(p1, p0, self.zlo, self.zhi)
            
        if p0zsliced:
            self.DistLamPpointPZ(p0z)
        if p1zsliced:
            self.DistLamPpointPZ(p1z)
        self.DistLamPedgePZF(p0z, p1z)

    def DistLamPtrianglePZO(self, p0, p1, p2):
        assert p0.z <= p1.z <= p2.z
        if p0.z <= self.zlo <= p2.z:
            self.DistLamPedgePZF(*PtriSlice(p0, p1, p2, self.zlo))
        if p0.z <= self.zhi <= p2.z:
            self.DistLamPedgePZF(*PtriSlice(p0, p1, p2, self.zhi))

    def DistLamPtrianglePZ(self, p0, p1, p2):
        if p0.z <= p1.z:
            if p0.z <= p2.z:
                if p1.z <= p2.z:
                    self.DistLamPtrianglePZO(p0, p1, p2)
                else:
                    self.DistLamPtrianglePZO(p0, p2, p1)
            else:
                self.DistLamPtrianglePZO(p2, p0, p1)
        elif p0.z <= p2.z:
            self.DistLamPtrianglePZO(p1, p0, p2)
        elif p1.z <= p2.z:
            self.DistLamPtrianglePZO(p1, p2, p0)
        else:
            self.DistLamPtrianglePZO(p2, p1, p0)

