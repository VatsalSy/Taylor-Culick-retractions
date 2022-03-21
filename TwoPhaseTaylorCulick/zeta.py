# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids
# Last Update: Dec 30 2020

import numpy as np
import os
import subprocess as sp
import sys

def readingXc(filename):
    fp = open(filename, "r")
    temp1 = fp.read()
    temp2 = temp1.split("\n")
    tTemp, XcTemp, YcTemp, VcTemp = [], [], [], []
    for n1 in range(len(temp2)):
        temp3 = temp2[n1].split(" ")
        if temp3 == ['']:
            pass
        else:
            tTemp.append(float(temp3[0]))
            XcTemp.append(float(temp3[1]))
            YcTemp.append(float(temp3[2]))
            VcTemp.append(float(temp3[3]))
    t = np.array(tTemp)
    Xc = np.array(XcTemp)
    Yc = np.array(YcTemp)
    Vc = np.array(VcTemp)
    return t, Xc, Yc, Vc

def writeDissipation(filename, Zo, Ro):
    exe = ["./getDissipation", filename, name, str(Zo), str(Ro), str(L0), str(N), str(Ohf), str(Ohs)]
    p = sp.Popen(exe, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = p.communicate()
    # temp1 = stderr.decode("utf-8")
    # temp2 = temp1.split("\n")
    # temp3 = temp2[0].split(" ")

# ----------------------------------------------------------------------------------------------------------------------


nGFS = 10000
ci = int(sys.argv[1])
L0 = 75
N = 2000
Ohf, Ohs = [0.0003, 1.0]
tmax = 101

folder = '%4.4d_ZetaData' % ci # output folder
if not os.path.isdir(folder):
    os.makedirs(folder)


tw, Zc, Rc, Vc = readingXc("%4d_X_Vjet.dat" % ci);

for ti in range(0, len(tw), 10):
    t, Zo, Ro, Vo = tw[ti], Zc[ti], Rc[ti], Vc[ti]
    if t > tmax:
        break
    place = "intermediate/snapshot-%5.4f" % t
    name = "%s/%5.4f.dat" % (folder, t)
    writeDissipation(place, Zo, Ro)
    print("Done t = %f at (R,Z) = (%f, %f)" % (t, Ro, Zo))
