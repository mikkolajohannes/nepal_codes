import numpy as np
import sys
import copy
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors
from matplotlib.cm import get_cmap
import netCDF4 as nc4
from wrf import (getvar, to_np, vertcross, CoordPair,get_cartopy, latlon_coords, xy_to_ll)

# Open the NetCDF file
#path = "./perpend_data/"
#path = "/home/local/mikkolaj/github/mikkolajohannes/nepal/valleys/puhti_perpenddata/east/curve2/east_c2_"
path = "/home/local/mikkolaj/github/mikkolajohannes/nepal/valleys/puhti_perpenddata/ncop/mid/ncopmid_"
#file = "2014-12-17_12:00:00.nc"
file = sys.argv[1] + ".nc"
data = path + file

date = file.split("_")[0]
hour = file.split("_")[1].split(":")[0]
minute = file.split("_")[1].split(":")[1]

ncfile = nc4.Dataset(data)
z = ncfile.variables["hgt"]
slope_wind = np.array(ncfile.variables["w_slo"])
valley_wind =  np.array(ncfile.variables["w_val"])
W_cross = np.array(ncfile.variables["w_vert"])
T_cross =  np.array(ncfile.variables["temp"])
len_sts = ncfile.variables["sts_len"][:]

slope_wind[slope_wind > 1e6] = np.NaN
valley_wind[valley_wind > 1e6] = np.NaN
W_cross[W_cross > 1e6] = np.NaN
T_cross[T_cross > 1e6] = np.NaN


#---------------------------------------------
#---------------------------------------------
#---------------------------------------------

fig = plt.figure(figsize=(12,9))
# ax_slope = fig.add_subplot(2,2,1)
# ax_valley = fig.add_subplot(2,2,2)
# ax_T = fig.add_subplot(2,2,4)
# ax_W = fig.add_subplot(2,2,3)
# ax_slope = fig.add_subplot(4,1,1)
# ax_valley = fig.add_subplot(4,1,3)
# ax_T = fig.add_subplot(4,1,4)
# ax_W = fig.add_subplot(4,1,2)

ax_valley = fig.add_subplot(2,1,1)
ax_T = fig.add_subplot(2,1,2)
#axes = [ax_slope, ax_valley, ax_T, ax_W]
axes = [ax_valley, ax_T]

#----------------------------------------------------
#----------------------------------------------------
#----------------------------------------------------
#----------------------------------------------------
# z1, z2 = 0, -1
# x1, x2 = 0, -1
#
# #east1
# z1, z2 = 2, 50
# x1, x2 = 0, -1
#
# #east2
# z1, z2 = 2, 60
# x1, x2 = 0, -1

# #mid3
# z1, z2 = 2, 60
# x1, x2 = 9, 23

# #ncopbot
# z1, z2 = 2, 60
# x1, x2 = 0, -1

# #ncopmid
z1, z2 = 2, 50
x1, x2 = 0, -1

#ncoptop
# z1, z2 = 3, 46
# x1, x2 = 15, -9


bounds_wind = np.arange(-40,42,2)
# contourf_slope = ax_slope.contourf(slope_wind[z1:z2,x1:x2], bounds_wind ,cmap=get_cmap("seismic"),extend="both")
# #contourf_slope = ax_slope.contourf(slope_wind, bounds_wind ,cmap=get_cmap("seismic"),extend="both")
# plt.colorbar(contourf_slope,ax=ax_slope)
#
# s_contour = ax_slope.contour(slope_wind[z1:z2,x1:x2], bounds_wind, colors="k", alpha=0.7)
# #s_contour = ax_slope.contour(slope_wind, bounds_wind, colors="k", alpha=0.7)
#
# #magic trick to manipulate inline labels
# fmt = {}
# strs = np.array((s_contour.levels[::2])).astype(int).astype(str)
#
# for l, s in zip(s_contour.levels[::2], strs):
#     fmt[l] = s
#
# ax_slope.clabel(s_contour,s_contour.levels[::2],fmt=fmt,inline=1)

#----------------------------------------------------
#----------------------------------------------------
#----------------------------------------------------
#----------------------------------------------------

#contourf_valley = ax_valley.contourf(valley_wind, bounds_wind ,cmap=get_cmap("seismic"),extend="both")
contourf_valley = ax_valley.contourf(valley_wind[z1:z2,x1:x2], bounds_wind ,cmap=get_cmap("seismic"),extend="both")
#plt.colorbar(contourf_valley,ax=ax_valley)

n = 2
#v_contour = ax_valley.contour(valley_wind, bounds_wind, colors="k", alpha=0.7)
v_contour = ax_valley.contour(valley_wind[z1:z2,x1:x2], bounds_wind, colors="k", alpha=0.7)
slope_wind_q = copy.copy(slope_wind[z1:z2,x1:x2])
slope_wind_q[::n] = 0
slope_wind_q[slope_wind_q == 0] = np.NaN
W_cross_q = copy.copy(W_cross[z1:z2,x1:x2])
W_cross_q[::n] = 0
W_cross_q[W_cross_q == 0] = np.NaN
ax_valley.quiver(slope_wind_q,W_cross_q, scale=5, scale_units='inches', pivot="mid")

#magic trick to manipulate inline labels
fmt = {}
strs = np.array((v_contour.levels[::2])).astype(int).astype(str)

for l, s in zip(v_contour.levels[::2], strs):
    fmt[l] = s

ax_valley.clabel(v_contour,v_contour.levels[::2],fmt=fmt,inline=1)


# Add a title
#ax_slope.set_title("Slope [m/s]", {"fontsize" : 12})
#ax_valley.set_title("Valley [m/s]",{"fontsize" : 12})


#----------------------------------------------------
#----------W AND T-----------------------------------
#----------------------------------------------------

# W_bounds = np.arange(-10,11,2)
# #W_bounds = [-10,-8,-6,-4,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,4,6,8,10]
# # W_contourf = ax_W.contourf(W_cross, W_bounds, cmap=get_cmap("seismic"),extend="both")
# # W_contourf = ax_W.contourf(W_cross, W_bounds, cmap=get_cmap("seismic"),extend="both")
# W_contourf = ax_W.contourf(W_cross[z1:z2,x1:x2], W_bounds, cmap=get_cmap("seismic"),extend="both")
# W_contourf = ax_W.contourf(W_cross[z1:z2,x1:x2], W_bounds, cmap=get_cmap("seismic"),extend="both")
# plt.colorbar(W_contourf,ticks=W_bounds,ax=ax_W)
#
# #v_contour = ax_valley.contour(valley_wind, bounds_wind, colors="k", alpha=0.7)
# w_contour = ax_W.contour(W_cross[z1:z2,x1:x2], W_bounds, colors="k", alpha=0.7)
#
# #magic trick to manipulate inline labels
# fmt = {}
# strs = np.array((w_contour.levels[::1])).astype(int).astype(str)
#
# for l, s in zip(w_contour.levels[::1], strs):
#     fmt[l] = s
#
# ax_W.clabel(w_contour,w_contour.levels[::1],fmt=fmt,inline=1)


#ticks=bounds[1:len(bounds)-1][::2]
cmap_T = mpl.colors.ListedColormap(['white','white'])
bounds_T = range(260,340)
norm_T = mpl.colors.BoundaryNorm(bounds_T,cmap_T.N)
sm_T = plt.cm.ScalarMappable(cmap=cmap_T, norm=norm_T)
sm_T.set_array([])
# T_contourf = ax_T.contourf(to_np(T_cross),cmap=cmap_T)
# T_contourf = ax_T.contourf(to_np(T_cross),cmap=cmap_T)
# T_contour = ax_T.contour(to_np(T_cross),bounds_T,colors="black")
T_contourf = ax_T.contourf(to_np(T_cross[z1:z2,x1:x2]),cmap=cmap_T)
T_contourf = ax_T.contourf(to_np(T_cross[z1:z2,x1:x2]),cmap=cmap_T)
T_contour = ax_T.contour(to_np(T_cross[z1:z2,x1:x2]),bounds_T,colors="black")

#magic trick to manipulate inline labels
fmt = {}
strs = np.array((T_contour.levels[::2])).astype(int).astype(str)

for l, s in zip(T_contour.levels[::2], strs):
    fmt[l] = s

ax_T.clabel(T_contour,T_contour.levels[::2],fmt=fmt,inline=1)

# Y AXIS LABELS

#tick labels to match thousands of meters
#error less than 25 meters
ton = np.arange(0,10001,500)
ton = ton[ton > to_np(z[z1])]
ton = ton[ton < to_np(z[z2])]

indx = 0
for zi in z:
    if (zi-ton[0]) > 0 and np.abs(zi-ton[0]) < 50:
        indx_low = indx
    if (zi-ton[-1]) < 0 and np.abs(zi-ton[-1]) < 50:
        indx_up = indx
        break
    indx += 1

vert_ticks = range(indx_low,indx_up+2)
vert_ticks = vert_ticks[::10]

for axi in axes:
    axi.set_yticks(vert_ticks)
    axi.set_yticklabels(ton,fontsize=10)
    axi.set_ylabel("Height [m]",fontsize=10)
    axi.set_facecolor("sienna")


# # print(to_np(z[vert_ticks[::20]]))
# # print(ton)
#
# vert_vals = to_np(z[z1:z2])
# v_ticks = np.arange(vert_vals.shape[0])
# vert_vals = np.around(vert_vals,decimals=0)
# v_ticks = np.around(v_ticks,decimals=0)
# # ax_slope.set_yticks(v_ticks[::10])
# # ax_slope.set_yticklabels(vert_vals[::10], fontsize=8)
# ax_valley.set_yticks(v_ticks[::10])
# ax_valley.set_yticklabels(vert_vals[::10], fontsize=8)
# ax_T.set_yticks(v_ticks[::10])
# ax_T.set_yticklabels(vert_vals[::10], fontsize=8)
# ax_W.set_yticks(v_ticks[::10])
# ax_W.set_yticklabels(vert_vals[::10], fontsize=8)
#
# ax_T.set_ylabel("Height [m]", fontsize=10)
# ax_W.set_ylabel("Height [m]", fontsize=10)

#----------------------------------------------------
#----------------------------------------------------
#----------------------------------------------------
#----------------------------------------------------

# Add a title
if int(hour) < 19:
	if minute == "30":
		nepal_time = str(int(hour)+6) + ":15"
	else:
		nepal_time = str(int(hour)+5) + ":45"
elif int(hour) > 19:
	if minute == "30":
		nepal_time = str(int(hour)-18) + ":15"
	else:
		nepal_time = str(int(hour)-19) + ":45"

if hour == "18" and minute == "30":
	nepal_time = "00:15"
time = date + "_" + hour + ":" + minute + ":00"
print(time)

# Add a title
#ax_T.set_title("Potential temperature  [K]", {"fontsize" : 12})
#ax_W.set_title("W [m/s]", {"fontsize" : 12})
title = "Cross section " + str(len_sts) + " km // " + time + " (" + nepal_time + ")"
fig.text(0.5,0.95,title,{"fontsize" : 12},horizontalalignment='center')

#ax_hgt.set_title(time)
fname = "/home/local/mikkolaj/github/mikkolajohannes/nepal/valleys/puhti_plots/ncopmid/circ_" + sys.argv[1].replace(":","_").replace("-","_") + ".png"
#fname = "testi_puhti.pdf"
plt.savefig(fname)
