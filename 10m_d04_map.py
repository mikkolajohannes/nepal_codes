from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.cm import get_cmap
import numpy as np
import os
import cartopy.crs as crs
from cartopy.feature import NaturalEarthFeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.patches as mpatches
import matplotlib.colors as colors
import netCDF4 as nc4

from wrf import (to_np, getvar, smooth2d, get_cartopy, cartopy_xlim,
                 cartopy_ylim, latlon_coords, interplevel,ALL_TIMES)


ncfile = nc4.Dataset("/home/mikkolaj/github/mikkolajohannes/nepal/wrfout_d04_2014-12-19_03:00:00")
LAT = getvar(ncfile, "XLAT")
LON = getvar(ncfile, "XLONG")

# Open the NetCDF file
data = "d04_variables.nc"
ncfile = nc4.Dataset(data)


# Get the sea level pressure
#hgt = getvar(wrflist, "HGT",timeidx=ALL_TIMES)
UU =  np.array(ncfile.variables["U"])
VV =  np.array(ncfile.variables["V"])
zz =  np.array(ncfile.variables["z"])
PBLH =  np.array(ncfile.variables["PBLH"])
hgt =  np.array(ncfile.variables["hgt"])

# dz = []
# for ii in range(0,zz.shape[1]-1):
#     for jj in range(0,zz.shape[2]-1):
#         dz.append(np.mean(zz[0,ii,jj,:])-hgt[ii,jj])
#
# print(np.amax(dz))
# print(np.amin(dz))

valleys_x, valleys_y = [], []
vals = ["west","ncop","mid","east"]
val_names = ["a) Gaurishankar","b) Khumbu","c) Makalu","d) Kanchanjunga"]

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

# print(np.amax(UU))
# print(np.amax(VV))

#avg daytime wind components
UU_daytime = np.empty((UU.shape[1]-1,UU.shape[2]-1,5))
VV_daytime = np.empty((VV.shape[1]-1,VV.shape[2]-1,5))

zi = 0

for dd in range(5):
    for xx in range(UU.shape[1]-1):
        for yy in range(UU.shape[2]-1):
            UU_daytime[xx,yy,dd] = np.mean(UU[zi,xx,yy,dd*12+5:dd*12+12])
            VV_daytime[xx,yy,dd] = np.mean(VV[zi,xx,yy,dd*12+5:dd*12+12])

# print(np.amax(UU_daytime))
# print(np.amax(VV_daytime))

#print(UU_daytime)
#exit()

wind_speed = np.sqrt(to_np(UU_daytime)**2 + to_np(VV_daytime)**2)

#print(*wind_speed)

#----wind dir
u_norm = to_np(UU_daytime) / np.sqrt(to_np(UU_daytime)**2 + to_np(VV_daytime)**2)
v_norm = to_np(VV_daytime) / np.sqrt(to_np(UU_daytime)**2 + to_np(VV_daytime)**2)
#---------------------------------------------------PLOT

fig = plt.figure(figsize=(15,12))
# Set the GeoAxes to the projection used by WRF

#times = [1,24,48,72,97]
#times = [1,24,48,72,97]

titles = ["17 Dec 2014 12-15 LT AVG","18 Dec 2014 12-15 LT AVG","19 Dec 2014 12-15 LT AVG","20 Dec 2014 12-15 LT AVG","21 Dec 2014 12-15 LT AVG"]
axes = []

for i in range(5):
#    axes.append(fig.add_subplot(2,3,i+1,projection=cart_proj))
    axes.append(fig.add_subplot(2,3,i+1))
    ax = axes[-1]
    #--------------------------------------------------------------------------------------
    #--------------------------------------------------------------------------------------
    #--------------------------------------------------------------------------------------
    levels_h = range(0,9000,1000)
    contours = plt.contour(range(hgt.shape[1]),range(hgt.shape[0]), hgt, levels_h, linewidths=1,colors="black")

#    windspd_levels = np.arange(0,26,1.0)
#    windspd_levels = [0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10,15,20,25,30]
    windspd_levels = [0,1,2,3,4,5,6,7,8,9,10,15,20,25,30]

    contourf = ax.contourf(range(hgt.shape[1]-1),range(hgt.shape[0]-1), wind_speed[:,:,i],
                            windspd_levels,
                            norm=colors.PowerNorm(gamma=0.45),
                            cmap=get_cmap("turbo"))#,
#
#     #tuulen suunta
    a = 10
    plt.quiver(range(hgt.shape[1]-1)[::a],range(hgt.shape[0]-1)[::a], u_norm[::a,::a,i], v_norm[::a,::a,i],
                    headlength=4,headaxislength=4,color="black",zorder=4)


    ax.set_xlim(20,230)
    ax.set_ylim(30,240)
    #--------------------------------------------------------------------------------------
    #--------------------------------------------------------------------------------------
    #--------------------------------------------------------------------------------------

    for val in range(0,4):
        if val==3:
            ax.plot(valleys_x[val][0:135],valleys_y[val][0:135],linestyle="--",linewidth=1,color="white")
        else:
            ax.plot(valleys_x[val],valleys_y[val],linestyle="--",color="white",linewidth=1)

    ax.set_title(titles[i])

    lat_ticks = [8,62,117,173,229,284]
    lat_ticklabels = ["26.5","27","27.5","28","28.5","29"]
    ax.set_yticks(lat_ticks)
    ax.set_yticklabels(lat_ticklabels)

    lon_ticks = [38,87,136,185,234,283]
    lon_ticklabels = ["86","86.5","87","87.5","88","88.5"]
    ax.set_xticks(lon_ticks)
    ax.set_xticklabels(lon_ticklabels)




    #plt.title("wind speed at 10 m")
cbar = plt.colorbar(contourf, ax=axes, shrink=.33, ticks=[0,1,2,3,4,5,6,7,8,9,10,15,20,25,30],orientation="horizontal",anchor=(1.0,4.0))
cbar.ax.set_xlabel("25m wind speed [m/s]")
cbar.ax.set_xticklabels([0,"",2,"",4,"",6,"",8,"",10,15,20,25,30])

plt.savefig("wind_d04.pdf",bbox_inches = 'tight')
    #plt.savefig("wind_d02.png")
