import zlib
import struct


class PNGslice:
    def __init__(self, fname, width, height):
        self.width = width
        self.height = height
        self.fout = open(fname,"wb")
        self.fout.write(b"\x89" + "PNG\r\n\x1A\n".encode('ascii'))

        colortype = 0     # true gray image (no palette)
        bitdepth = 8      # with one byte per pixel (0..255)
        compression = 0   # zlib (no choice here)
        filtertype = 0    # adaptive (each scanline seperately)
        interlaced = 0    # no
        block = struct.pack("!I4sIIBBBBB", 13, "IHDR".encode('ascii'), width, height,
                           bitdepth, colortype, compression, filtertype, interlaced)
        self.fout.write(block)
        self.fout.write(struct.pack("!I", zlib.crc32(block[4:])&0xFFFFFFFF))

        self.compressor = zlib.compressobj()
        self.lcompressed = []
        self.iy = 0
        
        self.blackline = b"\0"*width
        self.whiteline = b"\xFF"*width

    def AddRow(self, zls):
        assert self.iy < self.height
        self.lcompressed.append(self.compressor.compress(b"\0"))
        cline = self.whiteline
        piz = 0
        for iz in zls:
            assert iz <= self.width
            self.lcompressed.append(self.compressor.compress(cline[:iz-piz]))
            piz = iz
            cline = self.whiteline if cline == self.blackline else self.blackline
        self.lcompressed.append(self.compressor.compress(cline[:self.width-piz]))
        self.iy += 1
        
    def done(self):
        assert self.iy == self.height
        self.lcompressed.append(self.compressor.flush())
        self.fout.write(struct.pack("!I", sum(map(len, self.lcompressed))))
        IDAT = "IDAT".encode('ascii')
        crc = zlib.crc32(IDAT)
        self.fout.write(IDAT)
        for c in self.lcompressed:
            crc = zlib.crc32(c, crc)
            self.fout.write(c)
        self.fout.write(struct.pack("!I", crc&0xFFFFFFFF))

        block = "IEND".encode('ascii')
        bcrc = zlib.crc32(block)
        self.fout.write(struct.pack("!I4sI", 0, block, bcrc&0xFFFFFFFF))

        self.fout.close()
        
"""
ps = PNGslice("/home/goatchurch/geom3d/bitmapslices/xgxg.png", 1000, 1000)
for i in range(1000):
    ps.AddRow([100, 200+i//3])
ps.done()
"""
    
