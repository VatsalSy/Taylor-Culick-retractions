# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids

import numpy as np
import os
import sys

nGFS = 10000
ci = int(sys.argv[1])
Rhor = float(sys.argv[2])
Ohd = float(sys.argv[3])
Ohs = float(sys.argv[4])

name = "%d_getEnergy.dat" % ci

if os.path.exists(name):
    print("File %s found! New data will be appended to the file" % name)
for ti in range(nGFS):
    t = 1.00 * ti
    place = "intermediate/snapshot-%5.4f" % t
    if not os.path.exists(place):
        print("File %s not found!" % place)
    else:
        exe = "./getEnergyInerCap %s %s %s %s %s %s" % (place, name, Rhor, Ohd, Ohs)
        os.system(exe)
        print(("Done %d of %d" % (ti, nGFS)))
