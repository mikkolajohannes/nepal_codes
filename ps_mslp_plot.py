import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors
import wrf
import netCDF4 as nc4


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
QVAPOR = np.array(ncfile.variables["QVAPOR"])
MSLP = np.array(ncfile.variables["MSLP"])
MSLP=MSLP*100.0


P = P[:,0,:]
PB = PB[:,0,:]
T = P[:]+PB[:]-MSLP
T = np.flipud(T)

#
#
xticks = range(0,5*24*2,12)
xlabels = ["06","12","18","00",
            "06","12","18","00",
            "06","12","18","00",
            "06","12","18","00",
            "06","12","18","00",]

fig = plt.figure(figsize=(12,9))

axes_T = []
axes_dT = []

for i in range(1,6):
    axes_T.append(fig.add_subplot(2,5,i))
    axes_dT.append(fig.add_subplot(2,5,i+5))

#for axi in axes:

dT = np.zeros((T.shape))
for i in range(0,T.shape[0]):
    dT[i,:] = T[i,:]-np.mean(T[i,:])


titles = []
for i in range(17,22):
    titles.append(str(i) + "/12/2014")

i=0
for axi in axes_T:
    axi.set_yticks(np.flip(range(0,167))[::20])
    axi.set_yticklabels([])
#    axi.set_xlim(50000,100000)
    axi.set_ylim(0,167)
    # axi.set_xticks([70000,80000,90000,100000])
    # axi.set_xticklabels([700,800,900,1000])
    axi.grid(alpha=0.3)
    axi.set_title(titles[i])
    i+=1

axes_T[2].set_xlabel("p-mslp [hPa] vert. level 0",fontsize=14)
axes_dT[2].set_xlabel("$(p-mslp)-\\overline{(p-mslp)}$ [hPa] vert. level 0",fontsize=14)
axes_T[0].set_ylabel("Along valley grid points")
axes_dT[0].set_ylabel("Along valley grid points")



cmap = plt.get_cmap("Set1")

for i in range(0,5):#
    h=0
    for t in range(i*24*2,(i+1)*24*2)[::6]:
        axes_T[i].plot(T[:,t],range(0,167),'x',markersize=0.4,color=cmap(h))
        axes_dT[i].plot(dT[:,t],range(0,167),markersize=0.4,color=cmap(h),label=str(int(0.5*t)))
        h+=1

axes_T[0].set_yticklabels([0,20,40,60,80,100,120,140,160],fontsize=6)
axes_dT[0].set_yticklabels([0,20,40,60,80,100,120,140,160],fontsize=6)


for axi in axes_dT:
    axi.set_yticks(np.flip(range(0,167))[::20])
    axi.set_yticklabels([])
    axi.set_xlim(-500,500)
    axi.set_ylim(0,167)
    axi.grid(alpha=0.3)
    axi.set_xticks([-500,0,500])
    axi.set_xticklabels([-5,0,5])



lgnd = axes_dT[0].legend(bbox_to_anchor=(4.5, -0.2),ncol=8)
for i in range(0,8):
    lgnd.legendHandles[i]._legmarker.set_markersize(6)



fname = "/home/local/mikkolaj/github/mikkolajohannes/nepal/valleys/ncop/jan21/figures/ps_mslp_ncop.pdf"
plt.savefig(fname,bbox_inches = 'tight')
