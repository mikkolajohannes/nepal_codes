import numpy as np
import sys, os
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors
from matplotlib.cm import get_cmap
import netCDF4 as nc4
from wrf import (getvar, to_np, vertcross, CoordPair,get_cartopy, latlon_coords, xy_to_ll)
from scipy import interpolate


# os.system("rm slope_wind.pdf")
# os.system("rm swc_*")

ncfile = nc4.Dataset("/home/mikkolaj/github/mikkolajohannes/nepal/wrfout_d04_2014-12-19_03:00:00")
LAT = getvar(ncfile, "XLAT")
LON = getvar(ncfile, "XLONG")

lats = np.empty(LAT.shape[0])
for yi in range(LAT.shape[0]):
    lats[yi] = LAT[yi,0]

lons = np.empty(LON.shape[1])
for xi in range(LON.shape[1]):
    lons[xi] = LON[0,xi]

earth_radius = 6371e3
deg2distance = earth_radius*np.pi/180

#read valley centerline grid points into valleys_x and valleys_y
#0=west, 1=ncop, 2=mid, 3=east
valleys_x, valleys_y = [], []
vals = ["west","ncop","mid","east"]
val_names = ["a) Gaurishankar","b) Khumbu","c) Makalu","d) Kanchanjunga"]

titles = ["(a) Gaurishankar, west slope","(b) Gaurishankar, east slope","(c) Khumbu, west slope","(d) Khumbu, east slope",
            "(e) Makalu, west slope","(f) Makalu, east slope","(g) Kanchanjunga, west slope","(h) Kanchanjunga, east slope"]

for val in vals:
    valley_x = []
    valley_y = []
    infile = "../../" + val + "/valley_" + val + "_21.txt"
    with open(infile) as data_file:
        for line in data_file:
            parts = line.split() # split line into parts
            valley_x.append(int(parts[1]))
            valley_y.append(int(parts[0]))

    valleys_x.append(valley_x)
    valleys_y.append(valley_y)

#write the center points of each segment in each valley
#(given in nth number of grid point on the center line)
segs = np.empty((4,4)) #[valley,seg,t-b] = [west:east,seg_bot:seg_top,bot:top]


segs[0,:] = np.flip((95, 70, 45, 25)) #west
segs[1,:] = np.flip((100, 75, 50, 35)) #ncop
segs[2,:] = np.flip((155, 120, 90, 50)) #mid
segs[3,:] = np.flip((115, 90, 55, 20)) #east

slopes_y = [182, 172,
            158, 158,
            134, 134,
            113, 105,
            ####
            150, 150,
            135, 135,
            111, 111,
            93, 84,
            ####
            160, 160,
            131, 120,
            96, 103,
            64, 64,
            ####
            130, 130,
            93, 103,
            78, 68,
            62, 50]

slopes_x = [60, 72,
            50, 70,
            36, 56,
            37, 55,
            ####
            99, 115,
            97, 117,
            97, 117,
            91, 103,
            ####
            165, 185,
            155, 174,
            143, 172,
            148, 162,
            ####
            210, 220,
            195, 205,
            185, 203,
            170, 175]


data = "slope_variables.nc"

ncfile = nc4.Dataset(data)
U = np.array(ncfile.variables["U"])
V = np.array(ncfile.variables["V"])
U10 = np.array(ncfile.variables["U10"])
V10 = np.array(ncfile.variables["V10"])

print(U.shape,V.shape)
zi_hs = ["25","90","190","300","450","630","860","1050","1250","1430","1620","1900",np.NaN,np.NaN,np.NaN,np.NaN,np.NaN]

xticks = range(0,5*24*2,12)
xlabels = ["06","12","18","00",
            "06","12","18","00",
            "06","12","18","00",
            "06","12","18","00",
            "06","12","18","00",]
xlabels_date = ["06","12\n17/12/2014","18","00",
            "06","12\n18/12/2014","18","00",
            "06","12\n19/12/2014","18","00",
            "06","12\n20/12/2014","18","00",
            "06","12\n21/12/2014","18","00",]

z_max = 10

colors = ["blue","red","purple","orange"]

#for zi in range(0,12): #vertical levels
#for zi in range(-1,z_max+1): #vertical levels
for zi in [0]: #vertical levels
    title_counter = 0
    fig, axes = plt.subplots(4,2,figsize=(24,18))
    for vv in range(0,4): #valleys

        ax1 = axes[vv,0]
        ax2 = axes[vv,1]


        for ss in range(0,4):
            segi = segs[vv,ss] #segments

            slope_wind = np.empty((2,V.shape[2]))
            dx, dy = slopes_x[vv*4+ss*2+1]-slopes_x[vv*4+ss*2], slopes_y[vv*4+ss*2+1]-slopes_y[vv*4+ss*2]
            dx = np.cos(np.deg2rad(lats[slopes_y[vv*4+ss*2+1]]))*deg2distance*(lons[slopes_x[vv*4+ss*2+1]]-lons[slopes_x[vv*4+ss*2]])
            dy = deg2distance*(lats[slopes_y[vv*4+ss*2+1]]-lats[slopes_y[vv*4+ss*2]])

            if zi==-1:
                slope_wind[0,:] = (dx*U10[vv*8+ss*2,:]+dy*V10[vv*4+ss*2,:])/np.sqrt(dx**2+dy**2)
                slope_wind[1,:] = (dx*U10[vv*8+ss*2+1,:]+dy*V10[vv*4+ss*2+1,:])/np.sqrt(dx**2+dy**2)
            else:
                slope_wind[0,:] = (dx*U[vv*8+ss*2,zi,:]+dy*V[vv*4+ss*2,zi,:])/np.sqrt(dx**2+dy**2)
                slope_wind[1,:] = (dx*U[vv*8+ss*2+1,zi,:]+dy*V[vv*4+ss*2+1,zi,:])/np.sqrt(dx**2+dy**2)


            ax1.plot(range(0,V.shape[2]),-1.0*slope_wind[0,:],label=str(int(segi)),color=colors[ss])
            ax2.plot(range(0,V.shape[2]),slope_wind[1,:],label=str(int(segi)),color=colors[ss])


        ax1.set_title(titles[title_counter])
        title_counter += 1
        ax2.set_title(titles[title_counter])
        title_counter += 1

        for axi in [ax1,ax2]:
            for i in range(0,5):
                axi.fill_between([i*48+24,(i+1)*48],[12,12],[-12,-12],color="gray",alpha=0.2)

            axi.plot([0,240],[0,0],'k')
            axi.legend()
            axi.grid()
            axi.set_xticks(xticks)
            axi.set_xticklabels(xlabels_date)
            axi.set_yticks([-10,-7.5,-5,-2.5,0,2.5,5,7.5,10])
            axi.set_yticklabels(["-10","","-5","","0","","5","","10"])
            axi.set_xlim(24,240)
            axi.set_ylim(-10,10)

        ax1.set_ylabel("m s$^{-1}$",fontsize=14)
        ax2.set_ylabel("m s$^{-1}$",fontsize=14)

        if vv==3:
            props = dict(boxstyle='round', facecolor='white', alpha=.5, linewidth=2)
            ax1.text(0.5, -0.25, "Positive = up-slope wind\nNegative = down-slope wind", transform=ax1.transAxes, fontsize=18, verticalalignment='top',horizontalalignment="center", bbox=props)
            ax2.text(1.7, -0.25, "Positive = up-slope wind\nNegative = down-slope wind", transform=ax1.transAxes, fontsize=18, verticalalignment='top',horizontalalignment="center", bbox=props)

    if zi==-1:
        figname = "swc_10m.pdf"
    elif zi >-1 and zi < 10:
        figname = "swc_0" + str(zi) + ".pdf"
    else:
        figname = "swc_" + str(zi) + ".pdf"

#    figname = "pm5_test.pdf"
    plt.savefig(figname,bbox_inches = 'tight',dpi=300)

# unit_str = "pdfunite"
# for zi in range(0,z_max+1):
#     if zi < 10:
#         unit_str += " swc_0" + str(zi) + ".pdf"
#     else:
#         unit_str += " swc_" + str(zi) + ".pdf"
#
# unit_str += " slope_wind.pdf"
#
# os.system(unit_str)
