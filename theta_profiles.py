import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors
from matplotlib.cm import get_cmap
import netCDF4 as nc4

# Open the NetCDF file
#path = "./perpend_data/"
data = "/home/local/mikkolaj/github/mikkolajohannes/nepal/valleys/ncop/jan21/ncop_variables_50levels.nc"
#file = "2014-12-17_12:00:00.nc"

ncfile = nc4.Dataset(data)
z = np.array(ncfile.variables["z"])
T =  np.array(ncfile.variables["T"])
T = T+300
U = np.array(ncfile.variables["U"])
V = np.array(ncfile.variables["V"])
W = np.array(ncfile.variables["W"])
P = np.array(ncfile.variables["P"])
PB = np.array(ncfile.variables["PB"])
ALT = np.array(ncfile.variables["ALT"])
PBLH = np.array(ncfile.variables["PBLH"])



y = [20,40,60,80,100,140]
#y = [100]

bl_height = np.empty((len(y),PBLH.shape[1]))

for i in range(0,len(y)):
    bl_height[i,:] = z[y[i],0,:]+PBLH[y[i],:]+25

z_n = 10


fig = plt.figure(figsize=(24,4*len(y)))

axes = []
cmap = plt.get_cmap("Set1")

i=0
for yi in y:
    for d in range(0,5):
        c=0
        axes.append(fig.add_subplot(len(y),5,i*5+d+1))
        for h in range(0,24)[::4]:
            axes[-1].plot(T[yi,0:z_n,d*24*2+h*2],z[yi,0:z_n,d*24*2+h*2],color=cmap(c))
            c+=1

    for axi in axes[-5:-1]:
        axi.set_ylim(z[yi,0,10],z[yi,z_n,10])
    ylbl = "n=" + str(yi)
    axes[-5].set_ylabel(ylbl)
    i+=1

for axi in axes:
#    axi.set_xlim(290,330)
    axi.grid(color="gray",alpha=0.3)

# xticks = range(0,5*24*2,12)
#
# xlabels = ["06","12\n17/12/2014","18","00",
#             "06","12\n18/12/2014","18","00",
#             "06","12\n19/12/2014","18","00",
#             "06","12\n20/12/2014","18","00",
#             "06","12\n21/12/2014","18","00",]
#
#
#
# fig = plt.figure(figsize=(12,5*len(y)))
#
# axes = []
# for i in range(0,len(y)):
#     axes.append(fig.add_subplot(len(y),1,i+1))
#
#     # contourf = axes[-1].contourf(range(0,240),z[y[i],0:z_n,50],valley_massflux[i,0:z_n,:],range(-12,13),cmap=get_cmap("RdBu_r"))
#     # contour_0 = axes[-1].contour(range(0,240),z[y[i],0:z_n,50],valley_massflux[i,0:z_n,:],[-100.0,0.0,100.0],colors=["black","black","black"],alpha=0.3)
#     contourf = axes[-1].contourf(range(0,240),z[y[i],0:z_n,50],valley_wind[i,0:z_n,:],range(-15,16),cmap=get_cmap("RdBu_r"))
#     contour_0 = axes[-1].contour(range(0,240),z[y[i],0:z_n,50],valley_wind[i,0:z_n,:],[-100.0,0.0,100.0],colors=["black","black","black"],alpha=0.3)
#     axes[-1].plot(range(0,240),bl_height[i,:],":k",alpha=0.7)
#     fig.colorbar(contourf, ax=axes[-1])
#
#     txt = "n = " + str(y[i])
#     axes[-1].set_ylabel(txt)
#
#     axes[-1].set_xticks(xticks)
#     axes[-1].set_xticklabels(xlabels)
#     axes[-1].grid(color="grey",alpha=0.3,linestyle=":")
#
# axes[0].set_title("Massflux (positive=up-valley)  kg / (s m$^2$)")
# #contourf = ax.contourf(np.array((np.arange(0,z_n),np.arange(0,240))),z[y,0:z_n,:],V[y,0:z_n,:],range(-11,12),cmap=get_cmap("RdBu_r"))


fname = "./figures/theta_profiles.pdf"
plt.savefig(fname,bbox_inches = 'tight')
# fname = "./figures/massflux_timeseries.png"
# plt.savefig(fname,bbox_inches = 'tight')
