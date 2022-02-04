import math, re, sys
from collections import namedtuple


defaultfilenamein = "/home/goatchurch/geom3d/wheelcontour/AY.Cnc"
defaultfilenameout = "/home/goatchurch/geom3d/wheelcontour/AY-out.Cnc"

#defaultfilenamein = "c:\\pythong1g23\\ay.cnc"
#defaultfilenameout = "c:\\pythong1g23\\tyop.cnc"


class P2(namedtuple('P2', ['u', 'v'])):
    __slots__ = ()
    def __new__(self, u, v):
        return super(P2, self).__new__(self, float(u), float(v))
    def __repr__(self):
        return "P2(%s, %s)" % (self.u, self.v)
    def __add__(self, a):
        return P2(self.u + a.u, self.v + a.v)
    def __sub__(self, a):
        return P2(self.u - a.u, self.v - a.v)
    def __mul__(self, a):
        return P2(self.u*a, self.v*a)
    def __neg__(self):
        return P2(-self.u, -self.v)
    def __rmul__(self, a):
        raise TypeError
    def Lensq(self):
        return self.u*self.u + self.v*self.v
    def Len(self):
        if self.u == 0.0:  return abs(self.v)
        if self.v == 0.0:  return abs(self.u)
        return math.sqrt(self.u*self.u + self.v*self.v)
        
        
    def assertlen1(self):
        assert abs(self.Len() - 1.0) < 0.0001
        return True
        
    @staticmethod
    def Dot(a, b):
        return a.u*b.u + a.v*b.v

    @staticmethod
    def ZNorm(v):
        ln = v.Len()
        if ln == 0.0:  
            ln = 1.0
        return P2(v.u/ln, v.v/ln)
            
    @staticmethod
    def APerp(v):
        return P2(-v.v, v.u)
    @staticmethod
    def CPerp(v):
        return P2(v.v, -v.u)



try: 
    sendactivity
except:
    def sendactivity(*args, **kwargs):
        pass

sendactivity("clearalltriangles")
Dmaxarcrad = 200

def Main():
    fin, fout = GetFilestreams()
    gcodeheader, pts, gcodefooter = ParseFile(fin)
    sendactivity(contours=[pts], materialnumber=1)
    linarcseq = FindTanposFitCircles(pts, gap=15, endgap=3)
    gcodemiddle = GenerateGcodeMiddle(linarcseq)
    fout.write(gcodeheader)
    fout.write(gcodemiddle)
    fout.write(gcodefooter)
    if gcodefooter[-1] != '\n':
        fout.write("\n")
    fout.write(Getstats(pts, linarcseq))
    fout.write("\n")

def Getstats(pts, linarcseq):
    res = [ "; GCode file processed by curvefit2.py" ]
    rad0s = [ (linarc[0]-linarc[1]).Len()  for linarc in linarcseq  if len(linarc)==3 ]
    rad2s = [ (linarc[2]-linarc[1]).Len()  for linarc in linarcseq  if len(linarc)==3 ]
    raderr = min(abs(r0-r1)  for r0,r1 in zip(rad0s, rad2s))
    perr = [ ]
    nlines = sum(1  for linarc in linarcseq  if len(linarc)==2)
    narcs = sum(1  for linarc in linarcseq  if len(linarc)==3)
    res.append("; %d lines converted to %d lines and %d arcs" % (len(pts)-1, nlines, narcs))
    for p in pts:
        for linarc in linarcseq:
            p0, p2 = linarc[0], linarc[-1]
            if (p0.u<=p.u<=p2.u  if p0.u<=p2.u  else  p2.u<=p.u<=p0.u):
                break
        lam = (p.u - p0.u)/(p2.u - p0.u)
        assert -1e-5<=lam<=1+1e-5, (p, lam)
        if len(linarc) == 2:
            perr.append(abs(p0.v + lam*(p2.v - p0.v) - p.v))
        else:
            r = (p0 - linarc[1]).Len()
            perr.append(abs((p - linarc[1]).Len() - r))
    res.append("; to tolerance %.3F" % max(perr))
    return "\n".join(res)
    
    
def FindTanposFitCircles(pts, gap, endgap):
    tgens = [ ]
    tgens.append(tangency(pts[:endgap+1], -1))
    for i in range(0, len(pts), gap):
        tgens.append(tangency(pts[i:i+gap+1], 0))
    tgens.append(tangency(pts[-(endgap+1):], 1))
    sendactivity(contours=[[p,p+n]  for p, n in tgens])
    
    linarcgroups = [ ]
    for i in range(1, len(tgens)):
        (p0, n0), (p2, n2) = tgens[i-1], tgens[i]
        fbc = FitBiCircles(p0, n0, p2, n2)  # returns one segment or two arcs
        linarcgroups.append(((p0, n0), (p2, n2), fbc))
        
    linarcseq = [ ]
    for X0, X1, fbc in linarcgroups:
        linarcseq.extend(fbc)

    cconts = [ ]
    for linarc in linarcseq:
        if len(linarc) == 3:
            p0, c, p2 = linarc
            rad = (p0 - c).Len()
            rad2 = (p2 - c).Len()
            #assert abs(rad - rad2) < 1e-4, (rad, rad2)
            cont = [ ]
            for k in range(51):
                lam = k/50.0
                m = p0*(1-lam) + p2*lam
                #cont.append(m)
                cont.append(c + (m-c)*(rad/(m-c).Len()))
            #cconts.append([p1,p1+P2(0,1)])
            cconts.append(cont)
            
    sendactivity(contours=cconts, materialnumber=3)
    
    return linarcseq

# generate the G-code out and splice it back into the output file
def GenerateGcodeMiddle(linarcseq):
    outputlines = [ "G1 X%.3F Y%.3F" % linarcseq[0][0]]
    for linarc in linarcseq:
        if len(linarc) == 2:
            outputlines.append("G1 X%.3F Y%.3F" % linarc[1])
        else:
            assert len(linarc) == 3
            vc = linarc[1] - linarc[0]
            dc = 2 if P2.Dot(P2.CPerp(linarc[2]-linarc[0]), vc) > 0 else 3
            outputlines.append("G%d X%.3F Y%.3F I%.3F J%.3F" % (dc, linarc[2].u, linarc[2].v, vc[0], vc[1]))
    return '\n'.join(outputlines)
    
    
def ParseFile(fin):
    flist = [(line, re.search(".*?G1 X([\-\.e\d]+) Y([\-\.e\d]+)", line))  for line in fin ]
    g1iseqs = [ i  for i in range(len(flist)+1)  if (i == 0 or flist[i-1][1] is None) != (i == len(flist) or flist[i][1] is None) ]
    g1seqpairs = [(g1iseqs[j-1], g1iseqs[j])  for j in range(1, len(g1iseqs), 2)]
    assert g1seqpairs, "No sequence of G1 motions found in file"
    i0, i1 = max(g1seqpairs, key=lambda X:X[1]-X[0])
    pts = [ P2(float(m.group(1)), float(m.group(2)))  for l, m in flist[i0:i1] ]
    return "".join(line  for line, m in flist[:i0]), pts, "".join(line  for line, m in flist[i1:])


def tangency(pts, iend):
    # first apply linear regression driven by u-direction
    Nfac = 1.0/len(pts)
    su, sv = sum(p.u  for p in pts)*Nfac, sum(p.v  for p in pts)*Nfac
    su2, suv = sum(p.u**2  for p in pts)*Nfac, sum(p.u*p.v  for p in pts)*Nfac
    mb = (suv - su*sv) / (su2 - su**2)
    cb = (sv - mb*su)

    # use the error to the regression line to find the curvature and side
    errv = [(p.v - (mb*p.u+cb), p)  for p in pts]
    eev = (errv[0][0] + errv[-1][0])/2  # which way do the ends bend?
    if eev >= 0.0:
        tpt = min(errv, key=lambda X:X[0])[1]
        nvec = P2.ZNorm(P2(-mb,1))
    else:
        tpt = max(errv, key=lambda X:X[0])[1]
        nvec = P2.ZNorm(P2(mb,-1))
    if iend != 0:  # force tangent point to one of the ends
        tpt = pts[0] if iend == -1 else pts[-1]
    sendactivity(contours=[[tpt, tpt+nvec]])
    return tpt, nvec


DFitBiCircles = None
def FitBiCircles(p0, n0, p1, n1):
    global DFitBiCircles
    DFitBiCircles = (p0, n0, p1, n1)
    #(p0, n0, p1, n1) = DFitBiCircles
    
    p = P2(P2.Dot(n0, p1 - p0), P2.Dot(P2.APerp(n0), p1 - p0))
    n = P2(P2.Dot(n0, n1), P2.Dot(P2.APerp(n0), n1))
    # solve (p + n*s - (r,0))^2 = (r-s)^2
    # p^2 + s^2 + r^2 + 2p.n*s - 2p_u*r - 2n_u*r*s = r^2 - 2r*s + s^2
    # p^2 + 2p.n*s - 2p_v*r - 2n_v*r*s = -2r*s
    c, cs, cr, crs = p.Lensq(), 2*P2.Dot(p, n), -2*p.u, -2*n.u + 2
    # c + cs*s + cr*r + crs*r*s = 0
    # s = -(c + cr*r)/(cs + crs*r)

    # find at the halfway point where the 2 circles will be tangential
    nh = n + P2(1, 0)
    # solve vc = p+n*s - P2(r,0);  vc.APerp(nh) = 0
    nperp = P2.APerp(n)
    # 0 = (p+n*s - P2(r,0)) . (nperp + (0,1))
    #   = p.nperp - r*nperp_u + p_v + n_v*s
    #   = p.nperp + p_v - r*nperp_u + n_v*s
    pn = P2.Dot(p, nperp) + p.v
    #   = pn*(cs + crs*r) - r*nhperp_u*(cs + crs*r) - n_v*(c + cr*r)
    qa = -nperp.u*crs
    qb2 = (pn*crs - nperp.u*cs - n.v*cr)/2
    qc = pn*cs - n.v*c

    #pn*(cs + crs*r) - r*nperp.u*(cs + crs*r) - n.v*(c + cr*r)
    #s = -(c + cr*r)/(cs + crs*r)
    #P2.Dot(p,nperp) - r*nperp.u + p.v + n.v*s

    qdq = qb2**2 - qa*qc
    if qdq < 0 or qa == 0:
        return [(p0,p1)]
    qs = math.sqrt(qdq) / qa
    qm = -qb2 / qa
    q0 = qm + qs
    q1 = qm - qs
    qa*q1**2 + 2*qb2*q1 + qc

    qss = [ ]
    den0, den1 = (cs + crs*q0), (cs + crs*q1)
    if den0 != 0:
        qss.append((q0, -(c + cr*q0)/den0))
    if den1 != 0:
        qss.append((q1, -(c + cr*q1)/den1))
    qss = [(r,s)  for r,s in qss  if r>0 and s>0 and max(r,s) < Dmaxarcrad]
    qss.sort(key=lambda X:abs(X[0]-X[1]))
    for r, s in qss:
        c0, c1 = p0 + n0*r, p1 + n1*s
        assert abs((c0 - c1).Len() - abs(r - s)) < 1e-3 
        vc = P2.ZNorm(c1 - c0)
        m = c0 + vc*r*(+1 if P2.Dot(vc, n0)<0 else -1)
        assert abs((p0-c0).Len() - r) < 1e-4
        assert abs((p1-c1).Len() - s) < 1e-4
        rm, sm = (m-c0).Len(), (m-c1).Len()
        if abs(r - rm) > 1e-4 or abs(s - sm) > 1e-4:
            continue
        assert abs(r - rm) < 1e-4, (r, rm)
        assert abs(s - sm) < 1e-4, (s, sm)
        sendactivity(points=(c0, c1))
        sendactivity(contours=[(c0, c1)])
        sendactivity(points=[m], materialnumber=1)
        return [(p0, c0, m), (m, c1, p1)]
    return [(p0,p1)]
    

def GetFilestreams():
    if True or len(sys.argv) > 1:
        fnamein = "/home/goatchurch/geom3d/wheelcontour/AY.Cnc"
        fnameout = "/home/goatchurch/geom3d/wheelcontour/AY-out.Cnc"
        fin = open(fnamein)
        fout = open(fnameout, "w")
    else:
        fin = sys.stdin
        fout = sys.stout
    return fin, fout

Main()  # function call at end so we can work on code at top of file


