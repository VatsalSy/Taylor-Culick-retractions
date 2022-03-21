# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids
# Last Update: Dec 30 2020

import numpy as np
import os
import subprocess as sp
import matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib
from matplotlib.ticker import StrMethodFormatter
from matplotlib.ticker import FormatStrFormatter
from matplotlib.collections import LineCollection
import sys

matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['text.latex.preamble'] = [r'\boldmath']

def gettingFacets(filename, Tracer):
    if Tracer == 1:
        exe = ["./getFacet1", filename]
    else:
        exe = ["./getFacet2", filename]
    p = sp.Popen(exe, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = p.communicate()
    temp1 = stderr.decode("utf-8")
    temp2 = temp1.split("\n")
    segs = []
    skip = False
    if (len(temp2) > 1e2):
        for n1 in range(len(temp2)):
            temp3 = temp2[n1].split(" ")
            if temp3 == ['']:
                skip = False
                pass
            else:
                if not skip:
                    temp4 = temp2[n1+1].split(" ")
                    r1, z1 = np.array([float(temp3[1]), float(temp3[0])])
                    r2, z2 = np.array([float(temp4[1]), float(temp4[0])])
                    segs.append(((r1, z1),(r2,z2)))
                    skip = True
    return segs

def gettingTriplePoint(filename):
    exe = ["./getX0Y0V0", filename, name1, rOld, DistCutoff]
    p = sp.Popen(exe, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = p.communicate()
    temp1 = stderr.decode("utf-8")
    temp2 = temp1.split("\n")
    temp3 = temp2[0].split(" ")
    return float(temp3[0]), float(temp3[1]), float(temp3[2]), float(temp3[3])

def gettingXmYmVm(filename):
    exe = ["./getXmYmVm", filename, name2]
    p = sp.Popen(exe, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = p.communicate()
    temp1 = stderr.decode("utf-8")
    temp2 = temp1.split("\n")
    temp3 = temp2[0].split(" ")
    return float(temp3[0]), float(temp3[1]), float(temp3[2]), float(temp3[3])

# ----------------------------------------------------------------------------------------------------------------------


nGFS = 10000
ci = int(sys.argv[1])
Ldomain = int(sys.argv[2])

rOld, DistCutoff = sys.argv[3], sys.argv[4]
rminp, rmaxp, zminp, zmaxp = [Ldomain, 0, -Ldomain/4., Ldomain/4.]

name1 = "%4.4d_X0Y0V0.dat" % ci
name2 = "%4.4d_XmYmVm.dat" % ci

if os.path.exists(name1):
    print("File %s found! New data will be appended to the file" % name1)
if os.path.exists(name2):
    print("File %s found! New data will be appended to the file" % name2)
folder = 'TrackingTP' # output folder
if not os.path.isdir(folder):
    os.makedirs(folder)

for ti in range(nGFS):
    t = 100*ti
    place = "intermediate/snapshot-%5.4f" % t
    ImageName = "%s/%9.9d.png" %(folder, int(1e3*ti))
    if not os.path.exists(place):
        print("%s File not found!" % place)
    else:
        if os.path.exists(ImageName):
            print("%s Image present!" % ImageName)
        else:
            facets1 = gettingFacets(place, 1)
            facets2 = gettingFacets(place, 2)
            if (len(facets1) == 0 or len(facets2) == 0):
                print("Problem in the available file %s" % place)
            else:
                tp, zTP, rTP, vTP  = gettingTriplePoint(place)
                print("t %2.1f zTP %4.3f rTP %4.3f vTP %4.3e" % (tp, zTP, rTP, vTP))
                rOld = str(rTP)
                tp, zm, rm, vm  = gettingXmYmVm(place)
                print("t %2.1f zM %4.3f rM %4.3f vM %4.3e" % (tp, zm, rm, vm))
                ## Part to plot
                AxesLabel, TickLabel = [30, 25]
                fig, ax = plt.subplots()
                fig.set_size_inches(19.20, 10.80)
                rc('axes', linewidth=2)
                ## Drawing Facets
                line_segments1 = LineCollection(facets1, linewidths=2, colors='#fc8d59', linestyle='solid')
                ax.add_collection(line_segments1)
                line_segments2 = LineCollection(facets2, linewidths=2, colors='#1a9850', linestyle='solid')
                ax.add_collection(line_segments2)
                ax.plot([rTP], [zTP], 'bo')
                ax.plot([rm], [zm], 'ro')
                ax.plot([0, 0], [zminp, zmaxp],'--', color='grey')

                ax.set_xlabel(r'$\mathcal{R}$', fontsize=AxesLabel)
                ax.set_ylabel(r'$\mathcal{Z}$', fontsize=AxesLabel)
                ax.set_aspect('equal')
                ax.xaxis.set_major_formatter(FormatStrFormatter('$%.1f$'))
                ax.yaxis.set_major_formatter(FormatStrFormatter('$%.1f$'))
                ax.tick_params(labelsize=TickLabel)
                ax.set_xlim(rminp, rmaxp)
                ax.set_ylim(zminp, zmaxp)
                ax.set_title('$t = %4.3f$' % t, fontsize=AxesLabel)

                # plt.show()
                plt.savefig(ImageName,bbox_inches='tight')
                plt.close()
    print(("Done %d of %d" % (ti+1, nGFS)))
