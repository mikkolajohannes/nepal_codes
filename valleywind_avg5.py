import numpy as np
import sys, os
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors
from matplotlib.cm import get_cmap
import netCDF4 as nc4
from wrf import (getvar, to_np, vertcross, CoordPair,get_cartopy, latlon_coords, xy_to_ll)
from scipy import interpolate


os.system("rm valley_wind_avg.pdf")
os.system("rm vwc_*")

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
val_names = ["(a) Gaurishankar","(b) Khumbu","(c) Makalu","(d) Kanchanjunga"]

for val in vals:
    valley_x = []
    valley_y = []
    infile = "../../" + val + "/valley_" + val + "_21.txt"
    with open(infile) as data_file:
        for line in data_file:
            parts = line.split() # split line into parts
            valley_x.append(int(parts[1]))
            valley_y.append(int(parts[0]))

    # valleys_x.append(np.flip(valley_x))
    # valleys_y.append(np.flip(valley_y))
    valleys_x.append(valley_x)
    valleys_y.append(valley_y)

#write the center points of each segment in each valley
#(given in nth number of grid point on the center line)
segs = np.empty((4,4)) #[valley,seg,t-b] = [west:east,seg_bot:seg_top,bot:top]

# segs[0,:] = [80, 70, 45, 25] #west
# segs[1,:] = [100, 75, 50, 35] #ncop
# segs[2,:] = [145, 120, 90, 40] #mid
# segs[3,:] = [120, 90, 55, 20] #east

segs[0,:] = np.flip((95, 70, 45, 25)) #west
segs[1,:] = np.flip((100, 75, 50, 35)) #ncop
segs[2,:] = np.flip((155, 120, 90, 50)) #mid
segs[3,:] = np.flip((115, 90, 55, 20)) #east

# segs = np.empty((4,4,2)) #[valley,seg,t-b] = [west:east,seg_bot:seg_top,bot:top]
#
# segs[0,:,:] = [[90,105],[50,70],[30,50],[10,30]] #west
#
# segs[1,:,:] = [[95,120],[70,95],[50,70],[20,40]] #ncop
#
# segs[2,:,:] = [[140,160],[110,140],[60,100],[20,60]] #mid
#
# segs[3,:,:] = [[100,120],[70,90],[50,70],[5,25]] #east

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

z_max = 15

colors = ["blue","red","purple","orange"]

#for zi in range(0,12): #vertical levels
#for zi in range(0,z_max+1): #vertical levels
for zi in [3]: #vertical levels

    fig, axes = plt.subplots(4,1,figsize=(12,18))
    for vv in range(0,4): #valleys

        axi = axes[vv]

        # if vv==0:
        #     axi = axes[0,0]
        # elif vv==1:
        #     axi = axes[0,1]
        # elif vv==2:
        #     axi = axes[1,0]
        # elif vv==3:
        #     axi = axes[1,1]

    # fig, axes = plt.subplots(1,1,figsize=(12,5))
    # for vv in [0]: #valleys
    #     axi = axes

        axi.plot([0,240],[0,0],'k')
        for i in range(0,5):
            axi.fill_between([i*48+24,(i+1)*48],[20,20],[-25,-25],color="gray",alpha=0.2)

        axi.set_title(val_names[vv])

        data = "/home/local/mikkolaj/github/mikkolajohannes/nepal/valleys/2021_analysis/" + vals[vv] + "/" + vals[vv] + "_variables.nc"

        ncfile = nc4.Dataset(data)
        # z = np.array(ncfile.variables["z"])
        # T =  np.array(ncfile.variables["T"])
        # T = T+300
        U = np.array(ncfile.variables["U"])
        V = np.array(ncfile.variables["V"])
        # W = np.array(ncfile.variables["W"])
        # P = np.array(ncfile.variables["P"])
        # PB = np.array(ncfile.variables["PB"])
        # ALT = np.array(ncfile.variables["ALT"])
        # PBLH = np.array(ncfile.variables["PBLH"])
        for ss in range(0,4):
            segi = segs[vv,ss] #segments
            x, y = valleys_x[vv][int(segi)-4:int(segi)+5], valleys_y[vv][int(segi)-4:int(segi)+5]

            valley_wind = np.empty((5,V.shape[2]))
            t = 0
            for i in range(2,7):
#                dx_lon, dy_lat = x[i-2]-x[i+2], y[i-2]-y[i+2]
                dx = np.cos(np.deg2rad(lats[y[i]]))*deg2distance*(lons[x[i-2]]-lons[x[i+2]])
                dy = deg2distance*(lats[y[i-2]]-lats[y[i+2]])
                # print(lons[x[i-2]],lons[x[i+2]],lats[y[i-2]]-lats[y[i+2]])
                # print(dx,dy)
                valley_wind[:] = (dx*U[int(segi)-3+i,zi,:]+dy*V[int(segi)-3+i,zi,:])/np.sqrt(dx**2+dy**2)
                t+=1

            valley_wind_plot = np.empty((V.shape[2]))
            for t in range(V.shape[2]):
                valley_wind_plot[t] = np.mean(valley_wind[:,t])
            axi.plot(range(0,V.shape[2]),valley_wind_plot,label=str(int(segi)),color=colors[ss])

        axi.legend()
        axi.grid()
        axi.set_xticks(xticks)
        axi.set_xticklabels(xlabels_date)
#        axi.set_yticks([-5,-2.5,0,2.5,5,7.5,10])
        axi.set_yticks([-22.5,-20,-17.5,-15,-12.5,-10,-7.5,-5,-2.5,0,2.5,5,7.5,10])
#        axi.set_yticklabels(["-5","","0","","5","","10"])
        axi.set_yticklabels(["","-20","","-15","","-10","","-5","","0","","5","","10"])
        axi.set_ylabel("m s$^{-1}$",fontsize=14)
        axi.set_xlim(24,240)
        axi.set_ylim(-25,12.5)



    if zi < 10:
        figname = "vwc_0" + str(zi) + ".pdf"
    else:
        figname = "vwc_" + str(zi) + ".pdf"

#    figname = "pm5_test.pdf"

#    plt.savefig(figname,bbox_inches = 'tight')
    plt.savefig("test.pdf",bbox_inches = 'tight')

unit_str = "pdfunite"
for zi in range(0,z_max+1):
    if zi < 10:
        unit_str += " vwc_0" + str(zi) + ".pdf"
    else:
        unit_str += " vwc_" + str(zi) + ".pdf"

unit_str += " valley_wind_avg.pdf"

os.system(unit_str)




# valley_wind = np.empty((V.shape[0],V.shape[1],V.shape[2]))
# valley_massflux = np.empty((V.shape[0],V.shape[1],V.shape[2]))
# bl_height = np.empty((V.shape[0],PBLH.shape[1]))
#
# ddd = 2
# y_ll, y_ul = ddd,V.shape[0]-ddd
#for i in range(y_ll,y_ul):
#    dx, dy = valley_x[i+ddd]-valley_x[i-ddd], valley_y[i+ddd]-valley_y[i-ddd]
#    print(i,dx,dy)
#    dx, dy = 0, 1
#    valley_wind[i,:,:] = (dx*U[i,:,:]+dy*V[i,:,:])/np.sqrt(dx**2+dy**2)
#    valley_massflux[i,:,:] = valley_wind[i,:,:]/ALT[i,:,:]
#    bl_height[i,:] = z[i,0,:]+PBLH[i,:]+25
