from math import atan2, sin, cos, radians, ceil
from basicgeo import P3, P2, AlongAcc

class PostProcessor:
    def __init__(self, cfile, sendactivity):
        self.fout = open(cfile, "w")
        self.sendactivity = sendactivity
        self.zclear = 4
        self.fout.write("G21\n")
        self.fout.write("G0Z%.3f\n" % (self.zclear))
        self.x, self.y, self.z = 0, 0, self.zclear
        self.gcodemode = 0
        self.contcuts = [ ]
        self.contarcs = [ ]
        self.contlinks = [ ]

    def OutG(self, gcode, p, c):
        gline = [ ]

        if gcode == 2 or gcode == 3:
            gij = "I%.3fJ%.3f" % (c[0]-self.x, c[1]-self.y)
            if self.contcuts is not None:
                contarc = [ (self.x, self.y, self.z) ]
                p0 = P2(self.x, self.y)
                cc = P2(c[0], c[1])
                p1 = P2(p[0], p[1])
                v0 = p0 - cc
                v1 = p1 - cc
                v0perp = P2.CPerp(v0) if gcode == 2 else P2.APerp(v0)
                ang = atan2(P2.Dot(v1, v0perp), P2.Dot(v1, v0))
                ncc = int(ceil(ang/radians(5))+1)
                for i in range(1, ncc):
                    lam = i/ncc
                    lang = lam*ang
                    v = v0*cos(lang) + v0perp*sin(lang)
                    contarc.append(P3.ConvertGZ(cc + v, AlongAcc(lam, self.z, p[2])))
                contarc.append(p)
                if self.contarcs is not None:
                    self.contarcs.append(contarc)
                    self.contcuts.append([p])
                else:
                    self.contcuts[-1].extend(contarc[1:])
        else:
            gij = ""
            if self.contcuts is not None:
                self.contcuts[-1].append(p)
            
        if gcode != self.gcodemode:
            gline.append("G%d" % gcode)
            self.gcodemode = gcode
        if self.x != p[0]:
            gline.append("X%3f" % p[0])
            self.x = p[0]
        if self.y != p[1]:
            gline.append("Y%3f" % p[1])
            self.y = p[1]
        if self.z != p[2]:
            gline.append("Z%3f" % p[2])
            self.z = p[2]
        gline.append(gij)
        gline.append("\n")
        self.fout.write("".join(gline))

    def OutGseq(self, gseq, z):
        p0 = P3(gseq[0][0], gseq[0][1], z)
        self.OutRetractLink(p0)
        for gs in gseq[1:]:
            if gs[2] != 1 and max(abs(gs[0]-self.x), abs(gs[1]-self.y))>=0.05:
                self.OutG(gs[2], P3(gs[0], gs[1], z), P2(gs[3], gs[4]))
            else:
                self.OutG(1, P3(gs[0], gs[1], z), None)
        
    def OutRetractLink(self, p):
        lcont = [ (self.x, self.y, self.z) ]
        if self.z < self.zclear:
            self.fout.write("G0Z%.3f\n" % (self.zclear))
            lcont.append((self.x, self.y, self.zclear))
            self.fout.write("X%.3fY%.3f\n" % (p[0], p[1]))
            lcont.append((p[0], p[1], self.zclear))
        else:
            self.fout.write("G0X%.3fY%.3fZ%.3f\n" % (p[0], p[1], self.zclear))
            lcont.append((p[0], p[1], self.zclear))
        self.fout.write("G1Z%.3fF800\n" % p[2])
        self.x, self.y, self.z = p
        self.gcodemode = 1
        lcont.append((self.x, self.y, self.z))
        self.contlinks.append(lcont)
        self.contcuts.append([(self.x, self.y, self.z)])
            
    def OutLcont(self, cont):
        if cont[0] != (self.x, self.y, self.z):
            self.OutRetractLink(cont[0])
        for p in cont[1:]:
            self.fout.write("X%.3fY%.3fZ%.3f\n" % p)
        self.contcuts.append(cont)
        self.x, self.y, self.z = p
        

    def Close(self):
        self.contlinks.append([(self.x, self.y, self.z),(self.x, self.y, self.zclear)])
        self.fout.write("G0Z%.3f\n" % (self.zclear))
        self.fout.write("M2\n")
        self.fout.close()
        del self.fout
