import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors
from matplotlib.cm import get_cmap
import netCDF4 as nc4
from wrf import (getvar, to_np, vertcross, CoordPair,get_cartopy, latlon_coords, xy_to_ll)

# Open the NetCDF file
#path = "./perpend_data/"
path = "/home/local/mikkolaj/github/mikkolajohannes/nepal/valleys/ncop/jan21/along_theta/ncopalong_"
#file = "2014-12-17_12:00:00.nc"
file = sys.argv[1] + ".nc"
data = path + file

date = file.split("_")[0]
hour = file.split("_")[1].split(":")[0]
minute = file.split("_")[1].split(":")[1]

ncfile = nc4.Dataset(data)
z = ncfile.variables["hgt"]
T_cross =  np.array(ncfile.variables["temp"])
T_cross[T_cross > 1e6] = np.NaN

#---------------------------------------------
#---------------------------------------------
#---------------------------------------------

fig = plt.figure(figsize=(12,9))
ax_T = fig.add_subplot(1,1,1)

axes = [ax_T]

#----------------------------------------------------
#----------------------------------------------------
#----------------------------------------------------
#----------------------------------------------------
z1, z2 = 0, -2
x1, x2 = 0, -1
# z1, z2 = 0, -1
# x1, x2 = 0, -1

#ticks=bounds[1:len(bounds)-1][::2]

bounds_T = range(260,340)
whitewhite = [] #add "white" color for each theta contour
for i in bounds_T:
    whitewhite.append("white")
cmap_T = mpl.colors.ListedColormap(whitewhite)
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
ton = np.arange(0,10001,1000)
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
vert_ticks = vert_ticks[::20]

for axi in axes:
    axi.set_yticks(vert_ticks)
    axi.set_yticklabels(ton,fontsize=10)
    axi.set_ylabel("Height [m]",fontsize=12)
    axi.set_facecolor("sienna")

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
ax_T.set_title("Potential temperature  [K]", {"fontsize" : 12})
title = time + " (" + nepal_time + ")"
fig.text(0.5,0.95,title,{"fontsize" : 12},horizontalalignment='center')

#ax_hgt.set_title(time)
fname = "/home/local/mikkolaj/github/mikkolajohannes/nepal/valleys/ncop/jan21/figures/along/" + sys.argv[1].replace(":","_").replace("-","_") + ".pdf"
#fname = "testi_puhti.pdf"
plt.savefig(fname,bbox_inches = 'tight')
