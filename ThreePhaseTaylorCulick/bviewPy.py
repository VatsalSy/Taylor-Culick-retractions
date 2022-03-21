# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids

import numpy as np
import os
import sys

nGFS = 5000
Vmax = sys.argv[1]

folder = "Bview" # output folder
if not os.path.isdir(folder):
    os.makedirs(folder)
for ti in range(nGFS):
    t = 0.50 * ti
    place = "intermediate/snapshot-%5.4f" % t
    name = "%s/%8.8d.png" %(folder, int(1e3*t))
    if not os.path.exists(place):
        print("File %s not found!" % place)
    else:
        if os.path.exists(name):
            print("Image %s found!" % name)
        else:
            exe = "./getVideo %s %s %s" % (place, name, Vmax)
            os.system(exe)
            print(("Done %d of %d" % (ti, nGFS)))
