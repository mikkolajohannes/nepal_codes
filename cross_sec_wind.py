import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors
from matplotlib.cm import get_cmap
import cartopy.crs as crs
import scipy
import xarray as xr
import cartopy.feature as cfeature
from netCDF4 import Dataset
#from metpy.interpolate import cross_section
import wrf
from wrf import (getvar, get_basemap, to_np, vertcross, smooth2d, CoordPair, GeoBounds,
                 get_cartopy,interpline, latlon_coords, cartopy_xlim, cartopy_ylim)

from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

#python wrf_cross_sect.py wrfout_d04_2014-12-19_03:00:00

# Open the NetCDF file
#path = "/home/local/mikkolaj/codes/nepal/"
path = "/home/mikkolaj/github/mikkolajohannes/nepal/data/"
file = path + sys.argv[1]
ncfile = Dataset(file)

time = str(sys.argv[1])
time = time.replace("wrfout_d04_","")

# Get the WRF variables
z = getvar(ncfile, "z")
W = getvar(ncfile, "W")
hgt = getvar(ncfile,"HGT")
U = getvar(ncfile,"U")
V = getvar(ncfile,"V")

W = W[0:59,:,:]

valley_x = []
valley_y = []

with open("valley_ncop2.txt") as data_file:
    for line in data_file:
        parts = line.split() # split line into parts
        valley_x.append(int(parts[0]))
        valley_y.append(int(parts[1]))

start_point = CoordPair(x=valley_x[0], y=valley_y[0])
end_point = CoordPair(x=valley_x[-1], y=valley_y[-1])

#yksikkövektori
# unit_AB = np.empty((len(valley_x),2))
# for n in range(0,len(valley_x)-1):
#     AB_len = np.sqrt((valley_x[n+1]-valley_x[n])**2+(valley_y[n+1]-valley_y[n])**2)
#     AB = [valley_x[n+1]-valley_x[n],valley_y[n+1]-valley_y[n]]
#     unit_AB[n,0], unit_AB[n,1] = AB[0]/AB_len, AB[1]/AB_len

#print(unit_AB)
# apu = []
# U_nn = []
# V_nn = []
# for n in range(0,len(valley_x)-1):
#     dx = (valley_x[n+1]-valley_x[n])
#     dy = (valley_y[n+1]-valley_y[n])
#     U_n = U[:,valley_x[n],valley_y[n]]
#     V_n = V[:,valley_x[n],valley_y[n]]
#     U_nn.append(U[:,valley_x[n],valley_y[n]])
#     V_nn.append(V[:,valley_x[n],valley_y[n]])
#     apu.append((U_n*dx+V_n*dy)/(np.sqrt(U_n**2+V_n**2)*np.sqrt(dx**2+dy**2)))
#
# U_valley = U_nn*apu
# V_valley = V_nn*apu

valley_lat, valley_lon = wrf.xy_to_ll(ncfile,valley_y,valley_x)

valley = CoordPair(x=valley_x,y=valley_y)

interp_levels = range(0,int(hgt[valley_x[0],valley_y[0]])+500,50)

start_point = CoordPair(x=valley_y[0], y=valley_x[0])
end_point = CoordPair(x=valley_y[1], y=valley_x[1])


W_cross = vertcross(W,z,wrfin=ncfile,start_point=start_point,end_point=end_point,
                    levels=interp_levels,latlon=True)
U_cross = vertcross(U,z,wrfin=ncfile,start_point=start_point,end_point=end_point,
                    levels=interp_levels,latlon=True)
V_cross = vertcross(V,z,wrfin=ncfile,start_point=start_point,end_point=end_point,
                    levels=interp_levels,latlon=True)

for n in range(0,len(valley_x)-1,2) :
    start_point = CoordPair(y=valley_x[n], x=valley_y[n])
    end_point = CoordPair(y=valley_x[n+1], x=valley_y[n+1])
    W_cross_n = vertcross(W,z,wrfin=ncfile,
                        start_point=start_point,end_point=end_point,
                        levels=interp_levels, latlon=True)
    W_cross = xr.concat([W_cross,W_cross_n],dim="cross_line_idx")

    U_cross_n = vertcross(U,z,wrfin=ncfile,
                        start_point=start_point,end_point=end_point,
                        levels=interp_levels, latlon=True)
    U_cross = xr.concat([U_cross,U_cross_n],dim="cross_line_idx")

    V_cross_n = vertcross(V,z,wrfin=ncfile,
                        start_point=start_point,end_point=end_point,
                        levels=interp_levels, latlon=True)
    V_cross = xr.concat([V_cross,V_cross_n],dim="cross_line_idx")

valley3d = [range(0,59),range(0,len(valley_x)-1)]

# NCO-P station coordinates
ncop_y, ncop_x = 27.95, 86.82

# Get the lat/lon points
lats, lons = latlon_coords(hgt)

# Get the cartopy projection object
cart_proj = get_cartopy(hgt)

# Create a figure that will have 3 subplots
fig = plt.figure(figsize=(12,9))
#ax_hgt = fig.add_subplot(1,2,1,projection=cart_proj)
ax_W = fig.add_subplot(1,1,1)

#----------------------------------------------------
#----------------------------------------------------
#----------------------------------------------------

# levels = range(0,9000,2000)
#
# hgt_contours = ax_hgt.contour(to_np(lons), to_np(lats), to_np(hgt), levels, colors="black",
#             transform=crs.PlateCarree(),alpha=0.7)
#
# hgt_contourf = ax_hgt.contourf(to_np(lons), to_np(lats), to_np(hgt), levels,
#             transform=crs.PlateCarree(),
#             cmap=get_cmap("gist_earth"))
# cb_hgt = fig.colorbar(hgt_contourf, ax=ax_hgt, shrink=0.7,orientation="horizontal")
# cb_hgt.ax.tick_params(labelsize=8)
#
# hgt_contourf.set_clim(vmin=-4000,vmax=9000)
#
# gl = ax_hgt.gridlines(crs=crs.PlateCarree(), draw_labels=True,
#                   linewidth=2, color='gray', alpha=0.5, linestyle='--')
# gl.xlabels_top = False
# gl.ylabels_right = False
# gl.xformatter = LONGITUDE_FORMATTER
# gl.yformatter = LATITUDE_FORMATTER
#
# ax_hgt.plot(valley_lon,valley_lat,'yellow',linewidth=1.5,transform=crs.PlateCarree())
#
# ax_hgt.plot(ncop_x,ncop_y,'k*',markersize=12,transform=crs.PlateCarree())
# ax_hgt.quiver(to_np(valley_lon), to_np(valley_lat),
#            to_np(unit_AB[0,:]), to_np(unit_AB[:,0]),
#            transform=crs.PlateCarree(),headlength=2,headaxislength=2,color="black")

#----------------------------------------------------
#----------------------------------------------------
#----------------------------------------------------

W_contourf = ax_W.contourf(to_np(W_cross), cmap=get_cmap("jet"))
cb_W = fig.colorbar(W_contourf, ax=ax_W)

ax_W.contour(to_np(U_cross),cmap=get_cmap("Reds"))

# Set the x-ticks to use latitude and longitude labels
# coord_pairs = to_np(W_cross.coords["xy_loc"])
# x_ticks = np.arange(coord_pairs.shape[0])
# x_labels = [pair.latlon_str() for pair in to_np(coord_pairs)]
# ax_W.set_xticks([])
#
# # # Set the y-ticks to be height
vert_vals = to_np(W_cross.coords["vertical"])
v_ticks = np.arange(vert_vals.shape[0])
ax_W.set_yticks(v_ticks[::10])
ax_W.set_yticklabels(vert_vals[::10], fontsize=8)

# # Set the x-axis and  y-axis labels
ax_W.set_xlabel("Valley bottom path", fontsize=10)
ax_W.set_ylabel("Height (m)", fontsize=10)

# Add a title
ax_W.set_title("W, m/s", {"fontsize" : 10})

fname = "/home/mikkolaj/github/mikkolajohannes/nepal/valleys/figures/windcross_" + time + ".pdf"
plt.savefig(fname)
fname = "/home/mikkolaj/github/mikkolajohannes/nepal/valleys/figures/windcross_" + time + ".png"
plt.savefig(fname)
