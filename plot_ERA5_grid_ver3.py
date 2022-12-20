import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.cm import get_cmap
import numpy as np
import netCDF4 as nc4
import geopandas
import time
#import os
#import cartopy.crs as crs
#from cartopy.feature import NaturalEarthFeature
#from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.patches as mpatches
import cartopy.crs as ccrs
from cartopy.feature import NaturalEarthFeature
import cartopy.io.shapereader as shapereader
from mpl_toolkits.axes_grid1.inset_locator import inset_axes



data = nc4.Dataset("ERA5_d01_U.nc")
ERA_U = np.array(data.variables["u"]).astype(np.float)

data = nc4.Dataset("ERA5_d01_V.nc")
ERA_V = np.array(data.variables["v"]).astype(np.float)
lat, lon = np.array(data.variables["latitude"]).astype(np.float), np.array(data.variables["longitude"]).astype(np.float)

data = nc4.Dataset("ERA5_Nepal_17th21stDec2014_UV.nc")
dec_U = np.array(data.variables["u"]).astype(np.float)
dec_V = np.array(data.variables["v"]).astype(np.float)

#40 year avg
U400_avg = np.nanmean(ERA_U[1:40,10,:,:],0)
V400_avg = np.nanmean(ERA_V[1:40,10,:,:],0)

wspd400_avg = np.sqrt(U400_avg**2 + V400_avg**2)
U400_avg_norm = U400_avg / wspd400_avg
V400_avg_norm = V400_avg / wspd400_avg

wspd400 = np.empty((ERA_U.shape[0],ERA_U.shape[2],ERA_U.shape[3]))
for ti in range(10,ERA_U.shape[0]-10):
    wspd400[ti,:,:] = np.sqrt( ERA_U[ti,10,:,:]**2 + ERA_V[ti,10,:,:]**2 )

wspd400_std = np.std(wspd400,axis=0)
print(np.amin(wspd400_std),np.amax(wspd400_std))
#plot 40yr avg
fig = plt.figure(figsize=(15,8))
ax_6 = fig.add_subplot(2,3,1,projection=ccrs.PlateCarree())

contourf = ax_6.contourf(lon, lat, wspd400_avg, range(0,51,10), transform=ccrs.PlateCarree(),cmap=get_cmap("viridis"),zorder=1)

contour_std400 = ax_6.contour(lon, lat, wspd400_std, transform=ccrs.PlateCarree(),colors="black",zorder=3)
clabels = ax_6.clabel(contour_std400)


a = 5
ax_6.quiver(lon[::a], lat[::a], U400_avg_norm[::a,::a], V400_avg_norm[::a,::a],transform=ccrs.PlateCarree(),headlength=4,headaxislength=4,color="white",zorder=4)

ax_6.coastlines(color="gray",zorder=1)
country_borders = NaturalEarthFeature(
    category='cultural',
    name='admin_0_boundary_lines_land',
    scale='50m',
    facecolor='none')
ax_6.add_feature(country_borders, edgecolor='gray', linestyle="-", linewidth=1,zorder=1)

# get country borders
resolution = '10m'
category = 'cultural'
name = 'admin_0_countries'

shpfilename = shapereader.natural_earth(resolution, category, name)

# read the shapefile using geopandas
df = geopandas.read_file(shpfilename)

poly = df.loc[df['ADMIN'] == 'Nepal']['geometry'].values[0]

ax_6.add_geometries([poly], crs=ccrs.PlateCarree(), facecolor="none",edgecolor="r",zorder=2)

gl = ax_6.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
              linewidth=2, color='gray', alpha=0.1, linestyle='--')
gl.xlabels_top = False
gl.ylabels_right = False

ax_6.set_title("(a) Dec 1980-2019 ERA5 AVG")

titles = [" ","(b) 18 Dec 2014 ERA5 AVG","(c) 19 Dec 2014 ERA5 AVG","(d) 20 Dec 2014 ERA5 AVG","(e) 21 Dec 2014 ERA5 AVG"]
axes = []
iis = [0,1,2,4,5]
for i in range(1,5):

    #17-21 daily avg
    u_ti = np.nanmean(dec_U[i*4:i*4+4,:,:],0)
    v_ti = np.nanmean(dec_V[i*4:i*4+4,:,:],0)
    wspd_ti = np.sqrt(u_ti**2 + v_ti**2)
#    print(wspd_ti.shape)

    # diff_u_ti = u_ti - U400_avg
    # diff_v_ti = v_ti - V400_avg

    # u_stds = diff_u_ti / np.std(ERA_U[:,10,:,:],axis=0)

    #for wdir plot
    u_ti_norm = u_ti / wspd_ti
    v_ti_norm = v_ti / wspd_ti


    diff_wspd = wspd_ti - wspd400_avg #wspd difference to climatology
    how_many_stds = diff_wspd/wspd400_std #how many stand. dev. the difference is
    print(np.amin(how_many_stds),np.amax(how_many_stds))
    ax = fig.add_subplot(2,3,iis[i]+1,projection=ccrs.PlateCarree())
    ax.coastlines()

    ax.contourf(lon, lat, wspd_ti, range(0,51,10), transform=ccrs.PlateCarree(),cmap=get_cmap("viridis"),zorder=1)
#    contourf = ax.contourf(lon, lat, np.abs(diff_wspd), range(0,71,10), transform=ccrs.PlateCarree(),cmap=get_cmap("viridis"),zorder=1)

    ax.quiver(lon[::a], lat[::a], u_ti_norm[::a,::a], v_ti_norm[::a,::a],transform=ccrs.PlateCarree(),headlength=4,headaxislength=4,color="white",zorder=4)

    std_contour = ax.contour(lon, lat, how_many_stds, [-2,-1,0,1,2], transform=ccrs.PlateCarree(),colors="black",zorder=3)

    clabels = ax.clabel(std_contour)
    #
    # a = 5
    # plt.quiver(lons[::a,::a], lats[::a,::a],
    #            u_norm[times[i]+6,::a,::a], v_norm[times[i]+6,::a,::a],
    #            transform=crs.PlateCarree(),headlength=4,headaxislength=4,color="white",zorder=4)

    #--------------------------------------------------------------------------------------
    #--------------------------------------------------------------------------------------
    #--------------------------------------------------------------------------------------

    #plt.plot(ncop_x,ncop_y,'m*',markersize=12,transform=crs.PlateCarree())

    # Set the map bounds
    # ax.set_xlim(cartopy_xlim(hgt))
    # ax.set_ylim(cartopy_ylim(hgt))
    #
    # d04_box = ax.add_patch(mpatches.Rectangle(xy=[85.61807,26.434654], width=88.54175-85.61807, height=29.129082-26.434654,
    #                                 linewidth=2,
    #                                 facecolor='none',
    #                                 edgecolor="red",
    #                                 zorder=3,
    #                                 linestyle="-",
    #                                 transform=ccrs.PlateCarree()))
    #
    #
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=2, color='gray', alpha=0.1, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False

    ax.coastlines(color="gray",zorder=1)
    country_borders = NaturalEarthFeature(
        category='cultural',
        name='admin_0_boundary_lines_land',
        scale='50m',
        facecolor='none')
    ax.add_feature(country_borders, edgecolor='gray', linestyle="-", linewidth=1,zorder=1)
    ax.add_geometries([poly], crs=ccrs.PlateCarree(), facecolor="none",edgecolor="r",zorder=2)

    ax.set_title(titles[i])

    if i==3:
        axins = inset_axes(ax,
                            width="100%",
                            height="10%",
                            loc='lower center',
                            borderpad=-5
                           )
        fig.colorbar(contourf, cax=axins, orientation="horizontal",label="400 hPa wind speed [m s$^{-1}$]")


plt.savefig("fig_01_rev2.pdf",bbox_inches = 'tight',dpi=300)
