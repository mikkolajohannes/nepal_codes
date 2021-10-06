import numpy as np
import sys, os
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors
import netCDF4 as nc4
from matplotlib.lines import Line2D

#from mpltools import color

#os.system("rm *.pdf")

cmap = plt.get_cmap("Set1")

segs = np.empty((4,4)) #[valley,seg,t-b] = [west:east,seg_bot:seg_top,bot:top]

segs[0,:] = np.array((95, 70, 45, 25)) #west
segs[1,:] = np.array((100, 75, 50, 35)) #ncop
segs[2,:] = np.array((155, 120, 90, 50)) #mid
segs[3,:] = np.array((115, 90, 56, 20)) #east

colors = ["tab:red","tab:blue","tab:green","tab:purple","tab:orange","tab:olive","tab:brown","tab:pink"]
colors_crosses = ["blue","red","purple","orange"]

xticks = range(0,5*24*2,12)
xlabels = ["06","12","18","00",
            "06","12","18","00",
            "06","12","18","00",
            "06","12","18","00",
            "06","12","18","00",]

val = ["west","ncop","mid","east"]
val_names = ["a) Gaurishankar","b) Khumbu","c) Makalu","d) Kanchanjunga"]
#zis = [0,1,2,3,4,5,6,7,8,9,10,11]
zis = [0,4,8]

zi_hs = ["25","90","190","300","450","630","860","1050","1250","1430","1620","1900"]

fig = plt.figure(figsize=(12,15))
axes = []

plot_n = 0
for i in range(0,len(zis)):
    zi=zis[i]
    zi_h=zi_hs[i]


    for vali in range(4):
        valleyi = val[vali]

        axes.append(fig.add_subplot(3,4,vali+1+plot_n*4))
        data = "/home/local/mikkolaj/github/mikkolajohannes/nepal/valleys/2021_analysis/" + valleyi + "/" + valleyi + "_variables.nc"
        ncfile = nc4.Dataset(data)
        #z = np.array(ncfile.variables["z"])
        T =  np.array(ncfile.variables["T"])
        T = T+300
        T = T[:,zi,:]
        T = np.flipud(T)
        dT = np.zeros((T.shape[0],24*2))

        for yy in range(T.shape[0]):
            for tt in range(48):
                dT[yy,tt] = 0.5*(T[yy,3*24*2+tt]+T[yy,4*24*2+tt])-np.mean(T[yy,3*24*2:-1])

        h=0
        for t in range(0,24*2-1)[::6]:
            if vali==0:
                axes[-1].plot(dT[:,t],range(0,dT.shape[0]),markersize=0.4,color=colors[h],label=str(int(0.5*t)-24))
            else:
                axes[-1].plot(dT[:,t],range(0,dT.shape[0]),markersize=0.4,color=colors[h])

            h+=1

        if vali+plot_n*4 in [0,4,8]:
            axes[-1].set_ylabel("Vertical level "+str(zi),fontsize=14)


        if zi==0: axes[-1].set_title(val_names[vali],fontsize=14)



        #        xlab_T = "EAST // $\\theta$, [K] vert. level " + str(zi) + " (~" + str(zi_h) + " m)"
        #        xlab_dT = "EAST // $\\theta-\\overline{\\theta}$, [K] vert. level " + str(zi) + " (~" + str(zi_h) + " m)"
        # axes_T[2].set_xlabel(xlab_T ,fontsize=14)
        # axes_dT[2].set_xlabel(xlab_dT,fontsize=14)

        c=0
        for jj in range(0,4):
#            segi_x, segi_y = valley_nco_x[int(segs[i,jj])], valley_nco_y[int(segs[i,jj])]
#            segi_lat, segi_lon = wrf.xy_to_ll(ncfile, segi_y, segi_x)

            plt.scatter(4.75,segs[vali,jj],marker="x",s=45,color=colors_crosses[c])
            c+=1


    plot_n += 1


legend_elements = [Line2D([0], [0], color=colors[6], lw=2, label="00"),
                    Line2D([0], [0], color=colors[7], lw=2, label="03"),
                    Line2D([0], [0], color=colors[0], lw=2, label="06"),
                    Line2D([0], [0], color=colors[1], lw=2, label="09"),
                    Line2D([0], [0], color=colors[2], lw=2, label="12"),
                    Line2D([0], [0], color=colors[3], lw=2, label="15"),
                    Line2D([0], [0], color=colors[4], lw=2, label="18"),
                    Line2D([0], [0], color=colors[5], lw=2, label="21")]


lgnd = axes[0].legend(handles=legend_elements,bbox_to_anchor=(3.9, -2.55),ncol=8,title="Local time")


for axi in axes:
    axi.set_yticks(np.flip(range(0,190))[::20])
    axi.set_yticklabels([])
    axi.set_xlim(-5,5)
    axi.set_ylim(0,189)
    axi.grid(alpha=0.3)
    axi.set_xticks([-4,-2,0,2,4])
    axi.plot([0,0],[0,190],'k')

for j in [3,7,11]:
    axes[j].set_yticklabels(["","","Valley top","","","","","Valley\n entrance","","Plain"],rotation=270,fontsize=10)
    axes[j].yaxis.set_label_position("right")
    axes[j].yaxis.tick_right()



for j in [8,9,10,11]:
    axes[j].set_xlabel("$\\Delta $T [K]",fontsize=14)
    # for i in range(0,8):
    #     lgnd.legendHandles[i]._legmarker.set_markersize(6)

    # if zi < 10:
    #     fname = "T_0" + str(zi) + ".pdf"
    # else:
    #     fname = "T_" + str(zi) + ".pdf"

fname = "T_diurnal.pdf"

plt.savefig(fname,bbox_inches = 'tight')

#os.system("pdfunite T_*.pdf theta_diurnal.pdf")
