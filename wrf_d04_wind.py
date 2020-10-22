from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.cm import get_cmap
import numpy as np
import os,sys
import cartopy.crs as crs
from cartopy.feature import NaturalEarthFeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

from wrf import (to_np, getvar, smooth2d, get_cartopy, cartopy_xlim,
                 cartopy_ylim, latlon_coords)

# Open the NetCDF file
path = "../../../codes/nepal/"
file = path + sys.argv[1]
ncfile = Dataset(file)

time = str(sys.argv[1])
time = time.replace("wrfout_d04_","")

# Get the sea level pressure
hgt = getvar(ncfile, "HGT")
u10 = getvar(ncfile, "U10")
v10 = getvar(ncfile, "V10")

#-----itäinen laakso
# y_min = 144
# y_max = 232
# x_min = 44
# x_max = 135

#mt. everest vasen yläkulma
y_min = 50
y_max = 230
x_min = 30
x_max = 240

# y_min = 60
# y_max = 150
# x_min = 40
# x_max = 220

hgt = hgt[x_min:x_max,y_min:y_max]
u10 = u10[x_min:x_max,y_min:y_max]
v10 = v10[x_min:x_max,y_min:y_max]

# Get the latitude and longitude points
lats, lons = latlon_coords(hgt)

# NCO-P station coordinates
ncop_y, ncop_x = 27.95, 86.82

# Get the cartopy mapping object
cart_proj = get_cartopy(hgt)

#----wind speed
wind_speed = np.sqrt(u10**2 + v10**2)

# for i in range(0,np.abs(x_max-x_min)):
#     for j in range(0,np.abs(y_max-y_min)):
#         if hgt[i,j] > 6500:
#             wind_speed[i,j] = float("nan")
#             u10[i,j] = 0
#             v10[i,j] = 0

#----wind dir
u_norm = u10 / np.sqrt(u10**2 + v10**2)
v_norm = v10 / np.sqrt(u10**2 + v10**2)

#---------------------------------------------------PLOT

fig = plt.figure(figsize=(10,10))
# Set the GeoAxes to the projection used by WRF
ax = plt.axes(projection=cart_proj)
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
levels_h = range(0,8000,1000)
contours = plt.contour(to_np(lons), to_np(lats), to_np(hgt), levels_h,
            transform=crs.PlateCarree(),alpha=1,linewidths=1,colors="black")


cmap = mpl.colors.ListedColormap(['blue','dodgerblue','darkturquoise',
                               'cyan','greenyellow','yellow','orange',
                               'orangered','red','maroon','purple','grey'])
#cmap.set_over('grey')
#cmap.set_under('white')
bounds = [0,1,2,3,4,5,6,7,8,10,15,25,70]
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
contourf = plt.contourf(to_np(lons), to_np(lats), to_np(wind_speed),bounds,
            transform=crs.PlateCarree(),
            cmap=cmap,
            norm=norm)

plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),shrink=0.9,ax=ax,
                                   ticks=bounds[1:len(bounds)-1],extend='both')

#tuulen suunta
a = 5
plt.quiver(to_np(lons[::a,::a]), to_np(lats[::a,::a]),
           to_np(u_norm[::a,::a]), to_np(v_norm[::a,::a]),
           transform=crs.PlateCarree(),headlength=3,headaxislength=3,color="black")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------

plt.plot(ncop_x,ncop_y,'k*',markersize=12,transform=crs.PlateCarree())

# Set the map bounds
ax.set_xlim(cartopy_xlim(hgt))
ax.set_ylim(cartopy_ylim(hgt))


gl = ax.gridlines(crs=crs.PlateCarree(), draw_labels=True,
                  linewidth=2, color='gray', alpha=0.5, linestyle='--')
gl.xlabels_top = False
gl.ylabels_right = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER

title = "10m wind, " + time
plt.title(title)
fname = "/home/local/mikkolaj/github/mikkolajohannes/nepal/figures/wind_" + time + ".pdf"
plt.savefig(fname)
fname = "/home/local/mikkolaj/github/mikkolajohannes/nepal/figures/wind_" + time + ".png"
plt.savefig(fname)
