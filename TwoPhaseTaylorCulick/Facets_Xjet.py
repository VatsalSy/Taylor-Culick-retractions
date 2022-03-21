# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids

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

def gettingFacets(filename):
    exe = ["./getFacet", filename]
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
                    segs.append(((r1, z1),(r2, z2)))
                    segs.append(((r1, -z1),(r2, -z2)))
                    skip = True
    return segs

def gettingXjet(filename):
    exe = ["./getX_Vjet", filename, name]
    p = sp.Popen(exe, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = p.communicate()
    temp1 = stderr.decode("utf-8")
    temp2 = temp1.split("\n")
    temp3 = temp2[0].split(" ")
    return float(temp3[0]), float(temp3[1]), float(temp3[2]), float(temp3[3])

# ----------------------------------------------------------------------------------------------------------------------


nGFS = 4001
ci = int(sys.argv[1])
Ldomain = int(sys.argv[2])
tsnap = float(sys.argv[3])

rminp, rmaxp, zminp, zmaxp = [Ldomain, 0, -Ldomain/4., Ldomain/4.]

name = "%4.4d_X_Vjet.dat" % ci

if os.path.exists(name):
    print("File %s found! New data will be appended to the file" % name)
folder = 'TrackingJet' # output folder
if not os.path.isdir(folder):
    os.makedirs(folder)

for ti in range(nGFS):
    t = tsnap*ti
    place = "intermediate/snapshot-%5.4f" % t
    ImageName = "%s/%9.9d.png" %(folder, int(1e3*t))
    if not os.path.exists(place):
        print("%s File not found!" % place)
    else:
        if os.path.exists(ImageName):
            print("%s Image present!" % ImageName)
        else:
            facets = gettingFacets(place)
            if (len(facets) == 0):
                print("Problem in the available file %s" % place)
            else:
                tp, zTP, rTP, vTP  = gettingXjet(place)
                print("t %f zTP %7.6e rTP %7.6e vTP %7.6e" % (tp, zTP, rTP, vTP))
                ## Part to plot
                AxesLabel, TickLabel = [30, 25]
                fig, ax = plt.subplots()
                fig.set_size_inches(19.20, 10.80)
                rc('axes', linewidth=2)

                ## Drawing Facets
                line_segments = LineCollection(facets, linewidths=2, colors='#fc8d59', linestyle='solid')
                ax.add_collection(line_segments)

                ax.plot([rTP], [zTP], 'bo')
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
