from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.cm import get_cmap
import numpy as np
import os
import cartopy.crs as crs
from cartopy.feature import NaturalEarthFeature
import matplotlib.patches as mpatches
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from cartopy.feature import NaturalEarthFeature

import wrf
from wrf import (to_np, getvar, smooth2d, get_cartopy, cartopy_xlim,
                 cartopy_ylim, latlon_coords)

# Create a figure
def plot():
    fig = plt.figure(figsize=(20,8))
    # Set the GeoAxes to the projection used by WRF
    ax = fig.add_subplot(1,2,2,projection=cart_proj)
    #--------------------------------------------------------------------------------------
    #--------------------------------------------------------------------------------------
    #--------------------------------------------------------------------------------------
    levels = range(0,9000,1000)
    levels_500 = range(0,9000,200)

    # contours = plt.contour(to_np(lons), to_np(lats), to_np(hgt), levels, colors="black",
    #             transform=crs.PlateCarree(),alpha=0.4)
    # #plt.clabel(contours, inline=True, fontsize=8)
    # plt.contourf(to_np(lons), to_np(lats), to_np(hgt), 10,
    #             transform=crs.PlateCarree(),
    #             cmap=get_cmap("gist_earth"))
    # plt.clim(-6000,9000)
    # # Add a color bar
    # plt.colorbar(ax=ax, shrink=.98)


    vmin, vmax = -6000, 8000

    cmap = get_cmap("gist_earth")

    contours = ax.contour(to_np(lons), to_np(lats), to_np(hgt), levels, colors="black",
                transform=crs.PlateCarree(),alpha=0.4,linewidth=0.5,zorder=2)
    # contours_500 = plt.contour(to_np(lons), to_np(lats), to_np(hgt), levels_500, colors="black",
    #             transform=crs.PlateCarree(),alpha=0.4,linewidth=0.3,linestyle=":")

    #plt.clabel(contours, inline=True, fontsize=8)
    colorplot = ax.contourf(to_np(lons), to_np(lats), to_np(hgt), levels_500,
                transform=crs.PlateCarree(),vmin=vmin,vmax=vmax,
                cmap=cmap)
    ax.contourf(to_np(lons), to_np(lats), to_np(hgt), levels_500,
                transform=crs.PlateCarree(),vmin=vmin,vmax=vmax,
                cmap=cmap,zorder=1)
    # Add a color bar
    cbar = plt.colorbar(colorplot,ax=ax, shrink=.98)
    cbar.add_lines(contours)
    cbar.ax.set_ylabel("Surface height [m]")

    plt.plot(ncop_lon,ncop_lat,'k*',markersize=18,transform=crs.PlateCarree())

    gl = ax.gridlines(crs=crs.PlateCarree(), draw_labels=True,
                      linewidth=2, color='gray', alpha=0.5, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    #plt.show()



    for i in range(0,4):
        valley_nco_x = []
        valley_nco_y = []
        infile = "../../" + val[i] + "/valley_" + val[i] + "_21.txt"
        print("")
        print(infile)
        with open(infile) as data_file:
            for line in data_file:
                parts = line.split() # split line into parts
                valley_nco_x.append(int(parts[1]))
                valley_nco_y.append(int(parts[0]))
    #            print(to_np(LAT[valley_nco_y[-1],valley_nco_x[-1]]),to_np(LON[valley_nco_y[-1],valley_nco_x[-1]]))
    #    print("")
        # e20th_x, e20th_y = valley_nco_x[::20], valley_nco_y[::20]
        valley_nco_lat, valley_nco_lon = wrf.xy_to_ll(ncfile, valley_nco_x, valley_nco_y)
        # e20th_lat, e20th_lon = wrf.xy_to_ll(ncfile, e20th_y, e20th_x)


        plt.plot(valley_nco_lon,valley_nco_lat,color="yellow",transform=crs.PlateCarree(),zorder=3)
#        plt.scatter(e20th_lon,e20th_lat,color="yellow",transform=crs.PlateCarree())

        c=0
        for jj in range(0,4):
            segi_x, segi_y = valley_nco_x[int(segs[i,jj])], valley_nco_y[int(segs[i,jj])]
            segi_lat, segi_lon = wrf.xy_to_ll(ncfile, segi_x, segi_y)

            plt.scatter(segi_lon,segi_lat,marker="x",s=45,transform=crs.PlateCarree(),color=colors[c],zorder=4)
            c+=1

    # synop_lats = [28.63, 26.98, 27.3, 27.350, 26.48, 27.05, 27.07, 26.63]
    # synop_lons = [87.08, 87.35, 86.5, 87.670, 87.27, 88.26, 88.47, 88.32]
    #
    # for lati,loni in zip(synop_lats,synop_lons):
    #     plt.plot(loni,lati,'k+',markeredgewidth=2,markersize=18,transform=crs.PlateCarree())

    ncfile_d01 = Dataset("wrfout_d01_2014-12-17_00:00:00")

    # Get the sea level pressure
    hgt_d01 = getvar(ncfile_d01, "HGT")
    lats_d01, lons_d01 = latlon_coords(hgt_d01)


    ax2 = fig.add_subplot(1,2,1,projection=cart_proj)
    print(np.amin(wrf.to_np(hgt_d01)))
    hgt_d01 = np.array(hgt_d01)
    hgt_d01[hgt_d01==0] = np.NaN
    ax2.set_facecolor("skyblue")


    contours2 = ax2.contour(to_np(lons_d01), to_np(lats_d01), to_np(hgt_d01), levels, colors="black",
                transform=crs.PlateCarree(),alpha=0.4,linewidth=0.5)
    # contours_500 = plt.contour(to_np(lons), to_np(lats), to_np(hgt), levels_500, colors="black",
    #             transform=crs.PlateCarree(),alpha=0.4,linewidth=0.3,linestyle=":")

    #plt.clabel(contours, inline=True, fontsize=8)
    ax2.contourf(to_np(lons_d01), to_np(lats_d01), to_np(hgt_d01), levels_500,
                transform=crs.PlateCarree(),vmin=vmin,vmax=vmax,
                cmap=cmap)
    ax2.contourf(to_np(lons_d01), to_np(lats_d01), to_np(hgt_d01), levels_500,
                transform=crs.PlateCarree(),vmin=vmin,vmax=vmax,
                cmap=cmap)
    #ax2.clim(-6000,9000)


    gl = ax2.gridlines(crs=crs.PlateCarree(), draw_labels=True,
                      linewidth=2, color='gray', alpha=0.5, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

    d04_box = ax2.add_patch(mpatches.Rectangle(xy=[85.61807,26.434654], width=88.54175-85.61807, height=29.129082-26.434654,
                                    linewidth=2,
                                    facecolor='none',
                                    edgecolor="black",
                                    transform=crs.PlateCarree()))



    d03_box = ax2.add_patch(mpatches.Rectangle(xy=[80.0356,23.0023], width=92.749-80.0356, height=33.0705-23.0023,
                                    linewidth=2,
                                    facecolor='none',
                                    edgecolor="black",
                                    transform=crs.PlateCarree()))

    d02_box = ax2.add_patch(mpatches.Rectangle(xy=[72.4564,20.3042], width=95.1939-72.4564, height=38.1121-20.3042,
                                    linewidth=2,
                                    facecolor='none',
                                    edgecolor="black",
                                    transform=crs.PlateCarree()))

    ax2.text(85.8,26.6,"d04",transform=crs.PlateCarree(),fontsize="large")
    ax2.text(80.3,23.3,"d03",transform=crs.PlateCarree(),fontsize="large")
    ax2.text(72.5,20.5,"d02",transform=crs.PlateCarree(),fontsize="large")
    ax2.text(65.5,15.5,"d01",transform=crs.PlateCarree(),fontsize="large")

    ax.text(85.7,26.5,"d04",transform=crs.PlateCarree(),fontsize="x-large")

    ax2.plot(85.333,27.7,'r*',markersize=10,transform=crs.PlateCarree())
    ax2.plot(77.069710,28.679079,'r*',markersize=10,transform=crs.PlateCarree())

    #ax2.plot(ncop_lon,ncop_lat,'k*',markersize=18,transform=crs.PlateCarree())


    country_borders = NaturalEarthFeature(
        category='cultural',
        name='admin_0_boundary_lines_land',
        scale='50m',
        facecolor='none')
    ax2.add_feature(country_borders, edgecolor='k', linestyle="-", linewidth=1, alpha=0.8)
#    ax.add_feature(country_borders, edgecolor='k', linestyle="-", linewidth=1)


    ax2.set_title("(a)",fontsize=14)
    ax.set_title("(b)",fontsize=14)

    plt.savefig("valleys.pdf",bbox_inches = 'tight',dpi=300)

# Open the NetCDF file
ncfile = Dataset("/home/local/mikkolaj/github/mikkolajohannes/nepal/wrfout_d04_2014-12-19_09:00:00")

# Get the sea level pressure
hgt = getvar(ncfile, "HGT")
lap_hgt = getvar(ncfile, "LAP_HGT")
z = getvar(ncfile,"z")
LAT = getvar(ncfile,"XLAT")
LON = getvar(ncfile,"XLONG")

val = ["west","ncop","mid","east"]

# valley_nco_x = []
# valley_nco_y = []
# infile = "valley_mid_21.txt"
# with open(infile) as data_file:
#     for line in data_file:
#         parts = line.split() # split line into parts
#         valley_nco_x.append(int(parts[0]))
#         valley_nco_y.append(int(parts[1]))

# y_min, y_max = 20, 220
# x_min, x_max = 120, 200
#
#
# hgt = hgt[y_min:y_max,x_min:x_max]

# Get the latitude and longitude points
lats, lons = latlon_coords(hgt)

#print(lats,lons)

# NCO-P station coordinates
ncop_lat, ncop_lon = 27.95, 86.82

#synop observation coordinates


# Get the cartopy mapping object
cart_proj = get_cartopy(hgt)


# e20th_x, e20th_y = valley_nco_x[::20], valley_nco_y[::20]
#
# valley_nco_lat, valley_nco_lon = wrf.xy_to_ll(ncfile, valley_nco_y, valley_nco_x)
#
# e20th_lat, e20th_lon = wrf.xy_to_ll(ncfile, e20th_y, e20th_x)

# x0, y0 = 46, 125
# x1, y1 = 55, 144
# x2, y2 = 59, 167
# x3, y3 = 77, 186
# UV_x, UV_y = [y0,y0,y0,y1+3,y1,y1-3,y2,y2,y2,y3+3,y3,y3-3], [x0-5,x0,x0+4,x1-3,x1,x1+3,x2-5,x2,x2+5,x3-3,x3,x3+3]
# lat_UV_ncop, lon_UV_ncop = wrf.xy_to_ll(ncfile, UV_y, UV_x)

# for i in range(0,len(valley_nco_x)):
#     x,y = valley_nco_x[i], valley_nco_y[i]
#     pr = str(to_np(z[0,y,x])-to_np(hgt[y,x])) + " " + str(to_np(z[1,y,x])-to_np(hgt[y,x])) + " " + str(to_np(z[2,y,x])-to_np(hgt[y,x]))
#     print(pr)
#
# exit()

segs = np.empty((4,4)) #[valley,seg,t-b] = [west:east,seg_bot:seg_top,bot:top]

segs[0,:] = np.flip((95, 70, 45, 25)) #west
segs[1,:] = np.flip((100, 75, 50, 35)) #ncop
segs[2,:] = np.flip((155, 120, 90, 50)) #mid
segs[3,:] = np.flip((115, 90, 56, 20)) #east

colors = ["blue","red","purple","orange"]


plot()
