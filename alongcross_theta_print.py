import numpy as np
import sys, os
import xarray as xr
import netCDF4 as nc4
from wrf import (getvar, to_np, vertcross, CoordPair,get_cartopy, latlon_coords, xy_to_ll)

#from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

#python wrf_cross_sect.py wrfout_d04_2014-12-19_03:00:00

# Open the NetCDF file
path = "/home/local/mikkolaj/github/mikkolajohannes/nepal/"
file = path + sys.argv[1]
ncfile = nc4.Dataset(file)

time = str(sys.argv[1])
time = time.replace("wrfout_d04_","")

# Get the WRF variables
z = getvar(ncfile, "z")
hgt = getvar(ncfile, "HGT")
T = getvar(ncfile,"T")

T=T+300

valley_infile = "valley_ncop_jan21.txt"
valley_x = []
valley_y = []

with open(valley_infile) as data_file:
    for line in data_file:
        parts = line.split() # split line into parts
        valley_x.append(int(parts[0]))
        valley_y.append(int(parts[1]))

valley_lat, valley_lon = xy_to_ll(ncfile,valley_y,valley_x)

valley = CoordPair(x=valley_x,y=valley_y)

inter_min = float(0.95*hgt[valley_x[-1],valley_y[-1]])
inter_max = float(1.5*hgt[valley_x[0],valley_y[0]])
inter_step = 50

#---------------------------------------------------------------
#-----------CROSS SECTION---------------------------------------
#---------------------------------------------------------------

interp_levels = np.arange(inter_min,inter_max,inter_step)

start_point = CoordPair(x=valley_y[0], y=valley_x[0])
end_point = CoordPair(x=valley_y[1], y=valley_x[1])
T_cross = vertcross(T,z,wrfin=ncfile,start_point=start_point,end_point=end_point,
                    levels=interp_levels,latlon=True)

for n in range(1,len(valley_x)-1,2) :
    start_point = CoordPair(y=valley_x[n], x=valley_y[n])
    end_point = CoordPair(y=valley_x[n+1], x=valley_y[n+1])
    T_cross_n = vertcross(T,z,wrfin=ncfile,
                        start_point=start_point,end_point=end_point,levels=interp_levels, latlon=True)
    T_cross = xr.concat([T_cross,T_cross_n],dim="cross_line_idx")

outfile_name = "./along_theta/" + time + ".nc"
outfile = nc4.Dataset(outfile_name, 'w', format='NETCDF4')

outfile.createDimension("y", T_cross.shape[1]) #here x = from west slope to east slope
outfile.createDimension("z", T_cross.shape[0]) #interp_min and inter_max with 50m resolution
outfile.createDimension("longitude",1)
outfile.createDimension("latitude",1)

temp = outfile.createVariable("temp","f4",("z","y",),fill_value=np.NaN) #T_cross
hgt = outfile.createVariable("hgt","f4",("z",),fill_value=np.NaN) #height

hgt[:] = interp_levels      #m
temp[:] = T_cross           #K

outfile.close()
os.system("python alongcross_theta_plot.py " + time)
