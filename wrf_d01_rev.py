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

from wrf import (to_np, getvar, smooth2d, get_cartopy, cartopy_xlim,
                 cartopy_ylim, latlon_coords, interplevel,ALL_TIMES)

# Open the NetCDF file
wrflist = [Dataset("wrfout_d01_2014-12-17_00:00:00"), Dataset("wrfout_d01_2014-12-20_22:00:00")]



# Get the sea level pressure
hgt = getvar(wrflist, "HGT",timeidx=ALL_TIMES)
u = getvar(wrflist, "U",timeidx=ALL_TIMES)
v = getvar(wrflist, "V",timeidx=ALL_TIMES)
pb = getvar(wrflist, "PB",timeidx=ALL_TIMES)
p = getvar(wrflist, "P",timeidx=ALL_TIMES)

# print(u.shape)
# print(v.shape)
p_int = 40000.0

# print(p[100,10,100,100])
# print(pb[100,10,100,100])
# print(p[100,10,100,100]+pb[100,10,100,100])

u_int, v_int = interplevel(u[:,:,:,0:-1],p+pb,p_int), interplevel(v[:,:,0:-1,:],p+pb,p_int)

# Get the latitude and longitude points
lats, lons = latlon_coords(hgt)

# NCO-P station coordinates
ncop_y, ncop_x = 27.95, 86.82

# Get the cartopy mapping object
cart_proj = get_cartopy(hgt)

#----wind speed
wind_speed = np.sqrt(to_np(u_int)**2 + to_np(v_int)**2)


# for i in range(0,np.abs(x_max-x_min)):
#     for j in range(0,np.abs(y_max-y_min)):
#         if hgt[i,j] > 6500:
#             wind_speed[i,j] = float("nan")
#             u10[i,j] = 0
#             v10[i,j] = 0

#----wind dir
u_norm = to_np(u_int) / np.sqrt(to_np(u_int)**2 + to_np(v_int)**2)
v_norm = to_np(v_int) / np.sqrt(to_np(u_int)**2 + to_np(v_int)**2)


#---------------------------------------------------PLOT

fig = plt.figure(figsize=(10,12))
# Set the GeoAxes to the projection used by WRF

times = [1,24,48,72,97]
#times = [1,24,48,72,97]

titles = ["17 Dec 2014 12LT","(a) 18 Dec 2014 12LT","(b) 19 Dec 2014 12LT","(c) 20 Dec 2014 12LT","(d) 21 Dec 2014 12LT"]
axes = []

for i in range(1,5):
    axes.append(fig.add_subplot(2,2,i,projection=cart_proj))
    ax = axes[-1]
    #--------------------------------------------------------------------------------------
    #--------------------------------------------------------------------------------------
    #--------------------------------------------------------------------------------------
    levels_h = range(0,9000,1000)
    contours = plt.contour(to_np(lons), to_np(lats), to_np(hgt[times[i],:,:]), levels_h,
                transform=crs.PlateCarree(),linewidths=1,colors="black")


    #cmap = mpl.colors.ListedColormap(['navy','blue','dodgerblue','darkturquoise','cyan','greenyellow','yellow', 'orange','red'])
    #cmap.set_over('grey')
    #cmap.set_under('white')
    # bounds = [0.5,2,3,4,5,6,7,8,15,25]
    # norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    # contourf = plt.contourf(to_np(lons), to_np(lats), to_np(wind_speed),bounds,
    #             transform=crs.PlateCarree(),
    #             cmap=cmap,
    #             norm=norm)
    # plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),shrink=0.9,ax=ax,ticks=bounds,extend='both')

    contourf = ax.contourf(to_np(lons), to_np(lats), to_np(wind_speed[times[i]+6,:,:]),range(0,51,10),
                transform=crs.PlateCarree(),cmap=get_cmap("viridis"),zorder=1)

#    cbar = fig.colorbar(contourf,shrink=0.6)


    #tuulen suunta
    a = 5
    plt.quiver(to_np(lons[::a,::a]), to_np(lats[::a,::a]),
               to_np(u_norm[times[i]+6,::a,::a]), to_np(v_norm[times[i]+6,::a,::a]),
               transform=crs.PlateCarree(),headlength=4,headaxislength=4,color="white",zorder=4)

    #--------------------------------------------------------------------------------------
    #--------------------------------------------------------------------------------------
    #--------------------------------------------------------------------------------------

    #plt.plot(ncop_x,ncop_y,'m*',markersize=12,transform=crs.PlateCarree())

    # Set the map bounds
    ax.set_xlim(cartopy_xlim(hgt))
    ax.set_ylim(cartopy_ylim(hgt))

    d04_box = ax.add_patch(mpatches.Rectangle(xy=[85.61807,26.434654], width=88.54175-85.61807, height=29.129082-26.434654,
                                    linewidth=2,
                                    facecolor='none',
                                    edgecolor="red",
                                    zorder=3,
                                    linestyle="-",
                                    transform=crs.PlateCarree()))


    gl = ax.gridlines(crs=crs.PlateCarree(), draw_labels=True,
                      linewidth=2, color='gray', alpha=0.5, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

    ax.set_title(titles[i])


    #plt.title("wind speed at 10 m")
cbar = plt.colorbar(contourf, ax=axes, shrink=.5,orientation="horizontal",anchor=(0.5,1.7))
cbar.ax.set_xlabel("400 hPa wind speed [m s$^{-1}$]")

plt.savefig("test.pdf",bbox_inches = 'tight',dpi=300)
    #plt.savefig("wind_d02.png")
