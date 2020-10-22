from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.cm import get_cmap
import numpy as np
import os
import cartopy.crs as crs
from cartopy.feature import NaturalEarthFeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import wrf
from wrf import (to_np, getvar, smooth2d, get_cartopy, cartopy_xlim,
                 cartopy_ylim, latlon_coords)

# Create a figure
def plot():
    fig = plt.figure(figsize=(10,8))
    # Set the GeoAxes to the projection used by WRF
    ax = plt.axes(projection=cart_proj)
    #--------------------------------------------------------------------------------------
    #--------------------------------------------------------------------------------------
    #--------------------------------------------------------------------------------------
    levels = range(1000,8000,500)

    contours = plt.contour(to_np(lons), to_np(lats), to_np(hgt), levels, colors="black",
                transform=crs.PlateCarree(),alpha=0.4)
    #plt.clabel(contours, inline=True, fontsize=8)
    plt.contourf(to_np(lons), to_np(lats), to_np(hgt), 10,
                transform=crs.PlateCarree(),
                cmap=get_cmap("gist_earth"))
    plt.clim(-6000,9000)
    # Add a color bar
    plt.colorbar(ax=ax, shrink=.98)


    plt.plot(ncop_lon,ncop_lat,'k*',markersize=12,transform=crs.PlateCarree())

    plt.plot(coord_x0,coord_y0,'o',markersize=4,transform=crs.PlateCarree())
    plt.plot(valley_lon,valley_lat,'-',transform=crs.PlateCarree())
    # Set the map bounds
    ax.set_xlim(cartopy_xlim(hgt))
    ax.set_ylim(cartopy_ylim(hgt))

    # Add the gridlines
    #ax.gridlines(color="black", linestyle="dotted")

    plt.title("Title")

    gl = ax.gridlines(crs=crs.PlateCarree(), draw_labels=True,
                      linewidth=2, color='gray', alpha=0.5, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    #plt.show()
    plt.savefig("defvalley.pdf")
    plt.savefig("defvalley.png")


# Open the NetCDF file
ncfile = Dataset("/home/local/mikkolaj/codes/nepal/wrfout_d04_2014-12-19_12:00:00")

# Get the sea level pressure
hgt = getvar(ncfile, "HGT")
LAP = getvar(ncfile, "LAP_HGT")
LON = getvar(ncfile, "XLONG")
LAT = getvar(ncfile, "XLAT")

# Get the latitude and longitude points
lats, lons = latlon_coords(hgt)

# NCO-P station coordinates
ncop_lat, ncop_lon = 27.95, 86.82

# Get the cartopy mapping object
cart_proj = get_cartopy(hgt)

#starting point for valley definition

#keskim. laakso
y_start, x_start = 174, 145

#ncop
#y_start, x_start = 117,165
#y_start, x_start = 109,177
#vasen
#y_start, x_start = 60,170

#oikea laakso
#y_start, x_start = 215, 140
#y_start, x_start = 170, 56

"""

"""

coord_x0, coord_y0 = LON[x_start,y_start], LAT[x_start,y_start]

valley_y = []
valley_x = []

x = x_start
y = y_start

a = 1
for n in range(0,20):
    hgt_tmp = []
    for i in [-a,0,a]:
        for j in [-a,0,a]:
            if i==0 and j==0:
                continue
            if n > 4:
                if valley_x[-3] == x+i or valley_x[-2] == x+i:
                    if valley_y[-3] == y+j or valley_y[-2] == y+j:
                        continue
                else:
                    hgt_tmp.append(hgt[x+i,y+j])
                    if min(hgt_tmp) == hgt_tmp[-1]:
                        x_tmp, y_tmp = x+i, y+j
            else:
                hgt_tmp.append(hgt[x+i,y+j])
                if min(hgt_tmp) == hgt_tmp[-1]:
                    x_tmp, y_tmp = x+i, y+j

    valley_x.append(x_tmp)
    valley_y.append(y_tmp)
    x, y = x_tmp, y_tmp

# f  = open("valley_east.txt","a")
# for i in range(0,len(valley_x)-1):
#     f.write(str(valley_x[i]) + " " + str(valley_y[i]) + "\n")
# f.close()

#convert valley xy to coordinates
valley_lon = []
valley_lat = []




valley_lat, valley_lon = wrf.xy_to_ll(ncfile,valley_y,valley_x)

plot()
