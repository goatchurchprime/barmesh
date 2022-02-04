from PIL import Image, ImageDraw
import math

wi = 300
w1 = 200
ex = 4

sphs = [ {"x":0, "y":0, "z":0, "r":1, "c":0 }, 
         {"x":1.1, "y":0, "z":0.01, "r":0.3, "c":255 }
       ]

for zi in range(0, 100, 5):
    z = zi/50.0 - 1.0
    image = Image.new('L', (2*wi*ex, 2*wi*ex), color=255)
    draw = ImageDraw.Draw(image)
    for sph in sphs:
        cx, cy = sph["x"]*w1*ex + wi*ex, sph["y"]*w1*ex + wi*ex
        rsq = (sph["r"]*w1*ex)**2 - (z*w1*ex)**2
        if rsq > 0:
            r = math.sqrt(rsq)
            draw.ellipse((cx - r, cy - r, cx + r, cy + r), fill=sph["c"])

    simage = image.resize((2*wi, 2*wi), Image.ANTIALIAS)
    simage.save('/home/goatchurch/geom3d/sphereslices/slice%d.png' % zi)

