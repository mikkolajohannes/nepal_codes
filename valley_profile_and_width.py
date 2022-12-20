from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.cm import get_cmap
import numpy as np
import os
import wrf
from matplotlib.lines import Line2D




ncfile = Dataset("/home/mikkolaj/github/mikkolajohannes/nepal/wrfout_d04_2014-12-19_03:00:00")

# Get the sea level pressure
hgt = wrf.getvar(ncfile, "HGT")
lap_hgt = wrf.getvar(ncfile, "LAP_HGT")

valleys = ["west","ncop","mid","east"]
valleynames = ["(a) Gaurishankar","(b) Khumbu","(c) Makalu","(d) Kanchanjunga"]

segs = np.empty((4,4)) #[valley,seg,t-b] = [west:east,seg_bot:seg_top,bot:top]

segs[0,:] = np.flip((95, 70, 45, 25)) #west
segs[1,:] = np.flip((100, 75, 50, 35)) #ncop
segs[2,:] = np.flip((155, 120, 90, 50)) #mid
segs[3,:] = np.flip((115, 90, 56, 20)) #east

colors_crosses = ["blue","red","purple","orange"]


fig = plt.figure(figsize=(20,12))

for vali in range(len(valleys)):

    vals = ["Valley center","West ridge","East ridge"]
    path = "../../" + valleys[vali] + "/"
    infiles = [path + "valley_"+valleys[vali] + "_21.txt", path + "ridges/ridge_west.txt", path + "ridges/ridge_east.txt"]
    profiles = [[],[],[]]

    for i in range(3):

        valley_x = []
        valley_y = []
        infile = infiles[i]
        with open(infile) as data_file:
            for line in data_file:
                parts = line.split() # split line into parts
                if i==0:
                    valley_y.append(int(parts[0]))
                    valley_x.append(int(parts[1]))
                else:
                    valley_y.append(int(parts[1]))
                    valley_x.append(int(parts[0]))

        for y in range(len(valley_x)):
            profiles[i].append(wrf.to_np(float(hgt[valley_y[y],valley_x[y]])))

    profiles = np.array(profiles)

    width500 = []
    infile = "width_" + valleys[vali] + "_500.txt"
    with open(infile) as data_file:
        for line in data_file:
            parts = line.split() # split line into parts
            width500.append(float(parts[0]))


    width1000 = []
    infile = "width_" + valleys[vali] + "_1000.txt"
    with open(infile) as data_file:
        for line in data_file:
            parts = line.split() # split line into parts
            width1000.append(float(parts[0]))

    width1000 = np.array(width1000)
    width1000[width1000==0] = np.NaN

    width2000 = []
    infile = "width_" + valleys[vali] + "_2000.txt"
    with open(infile) as data_file:
        for line in data_file:
            parts = line.split() # split line into parts
            width2000.append(float(parts[0]))

    width2000 = np.array(width2000)

    widthBL = []
    infile = "width_" + valleys[vali] + "_BLmax.txt"
    with open(infile) as data_file:
        for line in data_file:
            parts = line.split() # split line into parts
            widthBL.append(float(parts[0]))

    widthBL = np.array(widthBL)


    ax = fig.add_subplot(2,2,vali+1)
    ax2 = ax.twinx()

    colors = ["k","r","r"]
    linestyles = ["-","--",":"]
    for i in range(1,3):
        ax.plot(profiles[i],color=colors[i],linestyle=linestyles[i],label=vals[i])

    depth_avg = []
    for yy in range(min(len(profiles[1]),len(profiles[2]))):
        depth_avg.append(0.5*(profiles[1][yy]+profiles[2][yy])-profiles[0][yy])

    depth_avg = np.array(depth_avg)
    ax2.plot(0.01*depth_avg,color="m",linestyle="-")
#    ax2.plot([0,180],[0.01*np.mean(depth_avg),0.01*np.mean(depth_avg)],"k:")
    ax2.plot([0,180],[0,0],color="gray",linestyle="-")

    ax.fill_between(range(len(profiles[0])),0,profiles[0],color="burlywood")

#    ax.set_ylim(0,8500)
    ax.set_ylim(0,14000)
#    ax.set_ylim(-6000,8500)

    ax.grid(linestyle=":",alpha=0.7)
    ax2.grid(linestyle=":",alpha=0.7)

    ax.set_title(valleynames[vali])
    ax.set_xlim(0,180)
    ax.set_xlabel("Along-valley grid points")
    ax.set_ylabel("Surface height [m]                                   ")

    #ax2.plot(range(len(profiles[0])),width500,label="width\n500m")
    ax2.plot(range(len(profiles[0])),width1000,color="b",linestyle="-",label="width\n1000m")
#    ax2.scatter(range(len(profiles[0])),width1000,color="b",linestyle="--",label="width\n1000m",alpha=0.5)

    #ax2.plot(range(len(profiles[0])),width2000,label="width\n2000m")
    #ax2.plot(range(len(profiles[0])),widthBL,label="width\nBLmax")


#    ax2.legend(loc=9)
    ax2.set_ylim(-70,50)
#    ax2.set_ylim(0,100)
    #ax2.set_ylabel("                                         Width [10km]\n                                         Depth [km]")
    ax2.text(186,20,"Width [x10 km]",color="blue",rotation=90,va="center")
    ax2.text(190,20,"Depth [km]",color="m",rotation=90,va="center")
    ax.set_yticks([0,2000,4000,6000])
    ax2.set_yticks([0,20,40])
    ax2.set_yticklabels([0,2,4])

#    ax2.plot(width1000)
    #ax2.plot(widthBL)
    #ax2.
    #ax2.set_ylim(0,150)

    c=0
    for jj in range(0,4):
        ax.scatter(segs[vali,jj],profiles[0][int(segs[vali,jj])]-50,marker="x",s=80,color=colors_crosses[c])
        c+=1



legend_elements = [Line2D([0], [0], color="burlywood", lw=6, label="Valley\ncenter"),
                    Line2D([0], [0], color=colors[1], linestyle=linestyles[1], lw=2, label="West\nridge"),
                    Line2D([0], [0], color=colors[2], linestyle=linestyles[2], lw=2, label="East\nridge"),
                    Line2D([0], [0], color="m", linestyle="-", lw=2, label="Depth"),
                    Line2D([0], [0], color="b", linestyle="-", lw=2, alpha=0.7, label="1000m\nwidth")]

ax.legend(handles=legend_elements,bbox_to_anchor=(0.3,-0.1),ncol=5,fontsize=12)

plt.savefig("valley_profiles.pdf",bbox_inches = 'tight')
