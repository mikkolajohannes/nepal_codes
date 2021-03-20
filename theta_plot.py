import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors
import netCDF4 as nc4
#from mpltools import color


# T = np.zeros((167,5*24*2)) #koordinaatit, aika
# j=0
# path = "/home/local/mikkolaj/github/mikkolajohannes/nepal/valleys/ncop/jan21/"
# H_infile = path + "jan21T_ncop.txt"
# with open(H_infile) as data_file:
#     for line in data_file:
#         parts = line.split() # split line into parts
#         for i in range(0,167):
#             T[i,j] = float(parts[i])+300
#         j=j+1


data = "/home/local/mikkolaj/github/mikkolajohannes/nepal/valleys/ncop/jan21/ncop_variables.nc"
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
PBLH = np.array(nc4.Dataset("ncop_PBLH.nc").variables["PBLH"])

T = T[:,0,:]


T = np.flipud(T)

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
    axi.set_xlim(285,325)
    axi.set_ylim(0,167)
    axi.set_xticks([290,300,310,320])
    axi.grid(alpha=0.3)
    axi.set_title(titles[i])
    i+=1

axes_T[2].set_xlabel("$\\theta$, [K] vert. level 0" ,fontsize=14)
axes_dT[2].set_xlabel("$\\theta - \overline{\\theta}$, [K] vert. level 0",fontsize=14)
axes_T[0].set_ylabel("Along valley grid points")
axes_dT[0].set_ylabel("Along valley grid points")

axes_T[0].set_yticklabels([0,20,40,60,80,100,120,140,160],fontsize=6)



cmap = plt.get_cmap("Set1")



for axi in axes_dT:
    axi.set_yticks(np.flip(range(0,167))[::20])
    axi.set_yticklabels([])
    axi.set_xlim(-5,5)
    axi.set_ylim(0,167)
    axi.grid(alpha=0.3)
    axi.set_xticks([-4,-2,0,2,4])
    axi.plot([0,0],[0,167],'k')

axes_dT[0].set_yticklabels([0,20,40,60,80,100,120,140,160],fontsize=6)
axes_dT[0].set_xlim(-10,5)
axes_dT[0].set_xticks([-10,-5,0,5])


for i in range(0,5):#
    h=0
    for t in range(i*24*2,(i+1)*24*2)[::6]:
        axes_T[i].plot(T[:,t],range(0,167),markersize=0.4,color=cmap(h))
        axes_dT[i].plot(dT[:,t],range(0,167),markersize=0.4,color=cmap(h),label=str(int(0.5*t)))
        h+=1

lgnd = axes_dT[0].legend(bbox_to_anchor=(4.5, -0.2),ncol=8)
for i in range(0,8):
    lgnd.legendHandles[i]._legmarker.set_markersize(6)






# #ax_1 = fig.add_subplot(3,1,1)
# # ax_2 = fig.add_subplot(3,1,2)
# # ax_3 = fig.add_subplot(3,1,3)
# #axes = [ax_1, ax_2, ax_3]
#
# axes = [ax_2,ax_3]
#
# for i in range(0,5):
#     for axi in axes:
#         axi.fill_between([i*48+24,(i+1)*48],[400,400],[-400,-400],color="gray",alpha=0.2)
#
#
# cmap = plt.get_cmap("RdPu_r")
#
# plots = []
#
# for i in range(0,167)[::20]:
#     ax_2.plot(T[i,:],color=cmap(i),label=str(i))
#     plots.append(ax_3.plot(T[i,:]-np.mean(T[i,:]),color=cmap(i),label=str(i)))
#
# # for i in range(20,167)[::20]:
# #     ax_2.plot(T[i,:],label=str(i))
# #     ax_3.plot(T[i,:]-np.mean(T[i,:]),label=str(i))
#
#
# ax_2.set_ylim(280,320)
# ax_3.set_ylim(-5,5)
# #-------------------
# #------------------1
# #-------------------
#
# #
# # plot_1west = ax_1.plot(H[6,:],'b')
# # plot_1mid = ax_1.plot(H[7,:],color="orange")
# # plot_1east = ax_1.plot(H[8,:],'g')
# #
# # #ax_1.set_ylim(300,320)
# #
# #
# # #-------------------
# # #------------------2
# # #-------------------
# #
# # plot_2west = ax_2.plot(H[3],'b')
# # plot_2mid = ax_2.plot(H[4],color="orange")
# # plot_2east = ax_2.plot(H[5],'g')
#
# #ax_2.set_ylim(297,307)
# #-------------------
# #------------------3
# #-------------------
#
# # ax_2.plot(H[1],'b')
# # ax_2.plot(H[4],color="orange")
# # ax_2.plot(H[7],'g')
# # ax_2.plot(Tp,'b--')
# # ax_2.set_ylim(285,310)
# # ax_2.set_ylabel("$\\theta$",fontsize=16)
# #
# #
# # bot = ax_3.plot(H[1]-np.mean(H[1]),'b')
# # mid = ax_3.plot(H[4]-np.mean(H[4]),color="orange")
# # top = ax_3.plot(H[7]-np.mean(H[7]),'g')
# # plain = ax_3.plot(Tp-np.mean(Tp),'b--')
# # #ax_3.set_ylim(293,303)
# # ax_3.set_xlabel("Local time",fontsize=16)
# # ax_3.set_ylim(-5,5.5)
# # ax_3.set_ylabel("$\\theta$-$\\theta_{mean}$",fontsize=16)
#
# for axi in axes:
#     axi.grid()
#     axi.set_xticks(xticks)
#     axi.set_xticklabels(xlabels)
#     axi.set_xlim(0,5*24*2)
#     # axi.set_ylim(-5,25)
#     axi.plot([0,240],[0,0],'k')
# #    axi.legend((plots[0],plots[int(0.5*len(plots))],plots[-1]),ncol=4)
#
# xlabels = ["06","12\n17/12/2020","18","00",
#             "06","12\n18/12/2020","18","00",
#             "06","12\n19/12/2020","18","00",
#             "06","12\n20/12/2020","18","00",
#             "06","12\n21/12/2020","18","00",]
# ax_3.set_xticklabels(xlabels)

# #ax_2.set_title("$\\theta$ [K]",fontsize=16)

# ax_3.legend((plots[160][0], plots[100][0],plots[60][0],plots[20][0]), ("Plain", 'Valley mouth', 'Middle valley', 'Valley head'),
#             loc='upper center', bbox_to_anchor=(0.5, -0.3),
#             fancybox=True, shadow=True, ncol=4,fontsize=14)


fname = "/home/local/mikkolaj/github/mikkolajohannes/nepal/valleys/ncop/jan21/figures/T_ncop.pdf"
plt.savefig(fname,bbox_inches = 'tight')
