# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids
# Last updated: 19-Nov-2020

import numpy as np
import os
import subprocess as sp
import matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib
from matplotlib.collections import LineCollection
from matplotlib.ticker import StrMethodFormatter
import sys

matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['text.latex.preamble'] = [r'\boldmath']

def readingXc(filename):
    fp = open(filename, "r")
    temp1 = fp.read()
    temp2 = temp1.split("\n")
    tTemp, XcTemp, YcTemp = [], [], []
    for n1 in range(1, len(temp2)):
        temp3 = temp2[n1].split(" ")
        if temp3 == ['']:
            pass
        else:
            tTemp.append(float(temp3[0]))
            XcTemp.append(float(temp3[1]))
            YcTemp.append(float(temp3[2]))
    t = np.array(tTemp)
    Xc = np.array(XcTemp)
    Yc = np.array(YcTemp)
    return t, Xc, -Yc

def gettingFacets(filename, tracer):
    print('Getting facets values')
    if tracer == 1:
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
                    segs.append(((r1, z1-0.005),(r2,z2-0.005)))
                    segs.append(((-r1, z1-0.005),(-r2,z2-0.005)))
                    skip = True
    return segs

# ----------------------------------------------------------------------------------------------------------------------


changeLeftAxes, Draft, tmaxOLD, XcmaxOLD = 1, int(sys.argv[5]), 1, 4
ci, Ldomain = int(sys.argv[1]), float(sys.argv[2])
tsnap, Rmin = float(sys.argv[3]), float(sys.argv[4])
tmin = tsnap
Zoom = 10
folder = 'TriplePointVideo'  # output folder

lw = 4
rWindow2, zWindow1, zWindow2 = [1, Ldomain/2, 1]
tw, Xc, Yc = readingXc("TriplePointTracking_%d_%3.2e.dat" % (ci, tsnap));
if changeLeftAxes:
    tmax, Rmax = max(tw), max(Xc)
else:
    tmax, Rmax = tmaxOLD, XcmaxOLD

if not os.path.isdir(folder):
    os.makedirs(folder)
for ti in range(len(tw)):
    t, Xo, Yo = tw[ti], Xc[ti], Yc[ti]
    if Draft:
        ti = len(tw)-1
        t, Xo, Yo = tw[ti], Xc[ti], Yc[ti]
    rmin1, rmax1, zmin1, zmax1 = [Ldomain, 0, -zWindow1/2, zWindow1/2]
    rmin2, rmax2, zmin2, zmax2 = [Xo+rWindow2/2, Xo-rWindow2/2, Yo-zWindow2/2, Yo+zWindow2/2]
    if rmin2 < 0:
        rmin2, rmax2 = [0, rWindow2]
    if rmax2 > Ldomain:
        rmin2, rmax2 = [Ldomain-rWindow2, Ldomain]
    place = "intermediate/snapshot-%5.4f" % t
    if Draft:
        name = "ci_%d.pdf" % ci
    else:
        name = "%s/%8.8d.png" %(folder, int(t*1000))

    if not os.path.exists(place):
        print("%s File not found!" % place)
    else:
        if os.path.exists(name) and not Draft:
            print("%s Image present!" % name)
        else:
            segs1 = gettingFacets(place,1)
            segs2 = gettingFacets(place,2)
            if (len(segs1) == 0 or len(segs2) == 0):
                print("Problem in the available file %s" % place)
            else:
                print("Doing %s" % place)
                # Part to plot
                AxesLabel, TickLabel = [40, 30]
                fig, (ax1, ax3) = plt.subplots(2,1, constrained_layout=True)
                fig.set_size_inches(19.20, 10.80)
                rc('axes', linewidth=2)
                facets1='#f46d43'
                facets2='#5e4fa2'

                ## Drawing Facets
                line_segments1 = LineCollection(segs1, linewidths=lw, colors=facets1, linestyle='solid')
                line_segments2 = LineCollection(segs2, linewidths=lw, colors=facets2, linestyle='solid')
                ax1.add_collection(line_segments2)
                ax1.add_collection(line_segments1)

                ax1.plot([Xo], [Yo], 'go')

                ax1.plot([rmin1, rmax1], [0.0, 0.0],'--',color='grey',linewidth=3)
                ax1.plot([rmin2, rmax2], [zmin2, zmin2],'-',color='grey',linewidth=1)
                ax1.plot([rmin2, rmax2], [zmax2, zmax2],'-',color='grey',linewidth=1)
                ax1.plot([rmin2, rmin2], [zmin2, zmax2],'-',color='grey',linewidth=1)
                ax1.plot([rmax2, rmax2], [zmin2, zmax2],'-',color='grey',linewidth=1)

                ax1.set_aspect('equal')
                ax1.set_xlim(rmin1, rmax1)
                ax1.set_ylim(zmin1, zmax1)
                ax1.tick_params(labelsize=TickLabel)
                ax1.set_xlabel(r'$\mathcal{R}$', fontsize=AxesLabel)
                ax1.set_ylabel(r'$\mathcal{Z}$', fontsize=AxesLabel)
                ax1.set_title('$t/t_\gamma = %4.3f$' % t, fontsize=TickLabel)


                ax1.yaxis.set_label_position("right")
                ax1.yaxis.tick_right()

                l, b, w, h = ax1.get_position().bounds
                left, bottom, width, height = [l+0.675*w, b+1.2*h-Zoom*zWindow2*h/Ldomain, Zoom*rWindow2*w/Ldomain, Zoom*zWindow2*h/Ldomain]
                ax2 = fig.add_axes([left, bottom, width, height])
                line_segments1 = LineCollection(segs1, linewidths=lw, colors=facets1, linestyle='solid')
                line_segments2 = LineCollection(segs2, linewidths=lw, colors=facets2, linestyle='solid')
                ax2.add_collection(line_segments2)
                ax2.add_collection(line_segments1)
                ax2.plot([Xo], [Yo], 'go')
                ax2.plot([rmin1, rmax1], [0.0, 0.0],'--',color='grey',linewidth=3)
                ax2.set_aspect('equal')
                ax2.set_xlim(rmin2, rmax2)
                ax2.set_ylim(zmin2, zmax2)
                labels = [item.get_text() for item in ax2.get_xticklabels()]
                empty_string_labels = ['']*len(labels)
                ax2.set_xticklabels(empty_string_labels)
                labels = [item.get_text() for item in ax2.get_yticklabels()]
                empty_string_labels = ['']*len(labels)
                ax2.set_yticklabels(empty_string_labels)

                ax3.plot(tw[1:ti+1], Xc[1:ti+1], color='#fdae61')
                tp, Xp = tw[1:ti+1], Xc[1:ti+1]
                ax3.plot(tp, Xp, color='#fdae61', marker='o', markersize=15, markeredgecolor='black',linestyle='none')

                # ax3.set_yscale('log')
                # ax3.set_xscale('log')

                ax3.tick_params(which='major', labelsize=TickLabel)
                ax3.tick_params(which='minor', labelsize=3*TickLabel/4)

                ax3.set_xlabel(r'$t/t_\gamma$', fontsize=AxesLabel)
                ax3.set_ylabel(r'$r$', fontsize=AxesLabel)

                ax3.plot([tmin, tmax], [Rmin, Rmin],'-',color='black',linewidth=2)
                ax3.plot([tmin, tmax], [Rmax, Rmax],'-',color='black',linewidth=2)
                ax3.plot([tmin, tmin], [Rmin, Rmax],'-',color='black',linewidth=2)
                ax3.plot([tmax, tmax], [Rmin, Rmax],'-',color='black',linewidth=2)

                ax3.set_xlim(tmin, tmax)
                ax3.set_ylim(Rmin, Rmax)
                ax3.set_aspect('equal')

                # plt.show()
                plt.savefig(name, bbox_inches="tight")
                plt.close()
                if Draft:
                    break
