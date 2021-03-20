import numpy as np
import sys, os
import wrf
import netCDF4 as nc4

var_name = sys.argv[1]

valley_infile = "/users/mikkolaj/nepal/valley_ncop_jan21.txt"
valley_x = []
valley_y = []

with open(valley_infile) as data_file:
    for line in data_file:
        parts = line.split() # split line into parts
        valley_x.append(int(parts[1]))
        valley_y.append(int(parts[0]))

DD =["17", "18", "19", "20", "21"]

hours = ["00", "01", "02", "03", "04", "05", "06",
        "07", "08", "09", "10", "11", "12",
        "13", "14", "15", "16", "17", "18",
        "19", "20", "21", "22", "23"]
mm = ["00", "30"]


outfile  = open("/scratch/project_2002790/Nepal_d04/massflux/ncop_variables.txt","w")

z_n = 30

U = np.empty([len(valley_x), z_n, len(DD)*len(hours)*len(mm)])
V = np.empty([len(valley_x), z_n, len(DD)*len(hours)*len(mm)])
W = np.empty([len(valley_x), z_n, len(DD)*len(hours)*len(mm)])
ALT = np.empty([len(valley_x), z_n, len(DD)*len(hours)*len(mm)])
T = np.empty([len(valley_x), z_n, len(DD)*len(hours)*len(mm)])
P = np.empty([len(valley_x), z_n, len(DD)*len(hours)*len(mm)])
PB = np.empty([len(valley_x), z_n, len(DD)*len(hours)*len(mm)])
z = np.empty([len(valley_x), z_n, len(DD)*len(hours)*len(mm)])

arrays = [U,V,W,ALT,T,P,PB,z]

var_name = ["U", "V", "W", "ALT", "T", "P", "PB", "z"]

t=0
for d in DD:
    for h in hours:
        for m in mm:
            file = "/scratch/project_2002790/Nepal_d04/wrfout_data/wrfout_d04_2014-12-" + d + "_" + h + ":" + m + ":00"
            ncfile = Dataset(file)

            for j in range(0,len(arrays))
                var = wrf.getvar(ncfile, var_name[j])

                for i in range(0,len(valley_x)):
                    x,y = valley_x[i], valley_y[i]

                    for zi in range(0,z_n):
                        arrays[j][i,zi,t] = wrf.to_np(var[zi,valley_y[i],valley_x[i]])
            t += 1


#outfile_name = "/scratch/project_2002790/Nepal_d04/massflux/ncop/" + time + ".nc"
outfile_name = "/scratch/project_2002790/Nepal_d04/massflux/ncop/" + time + ".nc"
outfile = nc4.Dataset(outfile_name, 'w', format='NETCDF4')

outfile.createDimension("y", len(valley_x))
outfile.createDimension("z", z_n)
outfile.createDimension("t", 5*24*2)
outfile.createDimension("longitude",1)
outfile.createDimension("latitude",1)

uwind = outfile.createVariable("U","f4",("y","z","t"),fill_value=np.NaN)
uwind[:] = U

vwind = outfile.createVariable("V","f4",("y","z","t"),fill_value=np.NaN)
vwind[:] = V

wwind = outfile.createVariable("W","f4",("y","z","t"),fill_value=np.NaN)
wwind[:] = W

dens = outfile.createVariable("ALT","f4",("y","z","t"),fill_value=np.NaN)
dens[:] = ALT

theta = outfile.createVariable("T","f4",("y","z","t"),fill_value=np.NaN)
theta[:] = T

base_pres = outfile.createVariable("PB","f4",("y","z","t"),fill_value=np.NaN)
base_pres[:] = PB

pres = outfile.createVariable("P","f4",("y","z","t"),fill_value=np.NaN)
pres[:] = P

hgt = outfile.createVariable("z","f4",("y","z","t"),fill_value=np.NaN) #height
hgt[:] = z

outfile.close()
