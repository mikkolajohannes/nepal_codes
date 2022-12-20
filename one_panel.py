import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors
from matplotlib.cm import get_cmap
import netCDF4 as nc4
import wrf


vals = ["west","ncop","mid","east"]
val_names = ["(a) Gaurishankar","(b) Khumbu","(c) Makalu","(d) Kanchanjunga"]

segs = np.empty((4,4)) #[valley,seg,t-b] = [west:east,seg_bot:seg_top,bot:top]

segs[0,:] = np.flip((95, 70, 45, 25)) #west
segs[1,:] = np.flip((100, 75, 50, 35)) #ncop
segs[2,:] = np.flip((155, 120, 90, 50)) #mid
segs[3,:] = np.flip((115, 90, 56, 20)) #east

ncfile = nc4.Dataset("/home/mikkolaj/github/mikkolajohannes/nepal/wrfout_d04_2014-12-19_03:00:00")
hgt = wrf.getvar(ncfile, "HGT")
LAT = wrf.getvar(ncfile, "XLAT")
LON = wrf.getvar(ncfile, "XLONG")

colors_crosses = ["blue","red","purple","orange"]

t1 = 3*24*2+6*2 #20th 6UTC -> 12LT
t2 = 3*24*2+9*2 #20th 9UTC -> 15LT
t3 = 3*24*2+3*2 #20th 6UTC -> 09LT
ts = [t1,t2,t3]
#ts = [t2]

for t in ts:
    fig, axes = plt.subplots(4,1,figsize=(10,4*5))

    for vv in range(0,4):
        val = vals[vv]
        val_name = val_names[vv]

        ax_1 = axes[vv]


        data = "/home/local/mikkolaj/github/mikkolajohannes/nepal/valleys/2021_analysis/"+val+"/"+val+"_variables.nc"

        ncfile = nc4.Dataset(data)
        z = np.array(ncfile.variables["z"])
        T =  np.array(ncfile.variables["T"])
        T = T+300
        U = np.array(ncfile.variables["U"])
        V = np.array(ncfile.variables["V"])
        # W = np.array(ncfile.variables["W"])
        # P = np.array(ncfile.variables["P"])
        # PB = np.array(ncfile.variables["PB"])
        # ALT = np.array(ncfile.variables["ALT"])
        PBLH = np.array(ncfile.variables["PBLH"])


        #######################################
        #valleywind
        #######################################


        valley_x = []
        valley_y = []
        infile = "../../"+val+"/valley_"+val+"_21.txt"
        with open(infile) as data_file:
            for line in data_file:
                parts = line.split() # split line into parts
                valley_x.append(int(parts[1]))
                valley_y.append(int(parts[0]))

        valley_lat = wrf.to_np(LAT[valley_y,0])
        valley_lon = wrf.to_np(LON[0,valley_x])

        # print(valley_lat)
        # print(valley_lon)

        ridges = [[],[]]

        for ii in [0,1]:

            ridge_x = []
            ridge_y = []
            if ii==0:
                infile = "../../"+val+"/ridges/ridge_west.txt"

            else:
                infile = "../../"+val+"/ridges/ridge_east.txt"

            with open(infile) as data_file:
                for line in data_file:
                    parts = line.split() # split line into parts
                    ridge_y.append(int(parts[1]))
                    ridge_x.append(int(parts[0]))

            for y in range(len(ridge_x)):
                ridges[ii].append(wrf.to_np(float(hgt[ridge_y[y],ridge_x[y]])))

        ridge_avg = []
        for yy in range(len(ridges[0])):
            ridge_avg.append(np.mean([ridges[0][yy],ridges[1][yy]]))

        valley_wind = np.empty((V.shape[0],V.shape[1],V.shape[2]))
        bl_height = np.empty((V.shape[0],PBLH.shape[1]))

        ddd = 5
        y_ll, y_ul = ddd,V.shape[0]-ddd

        for i in range(y_ll,y_ul):
            dx, dy = valley_x[i-ddd]-valley_x[i+ddd], valley_y[i-ddd]-valley_y[i+ddd]
            valley_wind[i,:,:] = (dx*U[i,:,:]+dy*V[i,:,:])/np.sqrt(dx**2+dy**2)
            bl_height[i,:] = z[i,0,:]+PBLH[i,:]+25

        valley_wind[0:y_ll,:,:],valley_wind[y_ul:-1,:,:] = valley_wind[y_ll,:,:], valley_wind[y_ul,:,:]


        z1, z2 = 0, 40
        y1, y2 = 0, -1

        z = z[:,z1:z2,:]

        y_plot = np.empty((z2,z.shape[0]))
        y_plot[:,:] = np.arange(0,V.shape[0])

#        wind_levels = [-15,-12,-10,-8,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,8,10,12,15]
        wind_levels = [-20,-15,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,15,20]
#        wind_levels = [-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10]

        avw_contourf = ax_1.contourf(y_plot, np.transpose(z[:,:,t]), np.transpose(valley_wind[:,z1:z2,t]),extend="both",
                                        norm=colors.SymLogNorm(linthresh=5.0, linscale=1.0,
                                              vmin=-15.0, vmax=15.0, base=10),
                                        levels=wind_levels,cmap=get_cmap("RdBu_r"))
    #    avw_contourf = ax_1.contourf(y_plot, np.transpose(z[:,:,t]), np.transpose(valley_wind[:,z1:z2,t]),wind_levels,cmap=get_cmap("RdBu_r"))
        contour_0 = ax_1.contour(y_plot, np.transpose(z[:,:,t]),np.transpose(valley_wind[:,z1:z2,t]),[-10.0,0.0,10.0],colors=["black","black","black"],alpha=0.3)
#        contour_15 = ax_1.contour(y_plot, np.transpose(z[:,:,t]),np.transpose(valley_wind[:,z1:z2,t]),[-15.0,15.0],colors=["black","black"],linestyle=":")


        #ax_1.set_title(val_name,fontsize=16)



        #######################################
        #theta
        #######################################

        T_contour = ax_1.contour(y_plot, np.transpose(z[:,:,t]), np.transpose(T[:,z1:z2,t]),np.arange(250,400,1.0),colors="black",alpha=0.6,linewidths=0.8)

        fmt = {}
        strs = np.array((T_contour.levels[::5])).astype(int).astype(str)

        for l, s in zip(T_contour.levels[::5], strs):
            fmt[l] = s

        ax_1.clabel(T_contour,T_contour.levels[::5],fmt=fmt,inline=1)

        #######################################
        #crosses
        #######################################

        path = "../../" + vals[vv] + "/"
        infile = path + "valley_"+ vals[vv] + "_21.txt"
        profiles = [[],[],[]]

        valley_x = []
        valley_y = []
        with open(infile) as data_file:
            for line in data_file:
                parts = line.split() # split line into parts
                valley_y.append(int(parts[0]))
                valley_x.append(int(parts[1]))

        for y in range(len(valley_x)):
            profiles[0].append(wrf.to_np(float(hgt[valley_y[y],valley_x[y]])))

        profiles = np.array(profiles)

        #######################################
        #rest
        #######################################


        if vv==3:
            cbar_ax = fig.add_axes([0.92, 0.32, 0.02, 0.35])
            cbar = fig.colorbar(avw_contourf,cax=cbar_ax,ticks=wind_levels,orientation="vertical")
            cbar.ax.set_yticklabels([-20,"",-10,"",-8,"",-6,"",-4,"",-2,"",0,"",2,"",4,"",6,"",8,"",10,"",20],fontsize=16)
            cbar.ax.set_ylabel("Along-valley wind [m s$^{-1}$]",fontsize=16)
            #cbar.add_lines(contour_0)

            #cbar = plt.colorbar(avw_contourf, ax=ax_1, shrink=.33, ticks=wind_levels,orientation="horizontal",anchor=(-1.0,-1.0))
#            cbar.ax.set_xlabel("25m wind speed [m/s]")
#            cbar.ax.set_xticklabels([0,"",2,"",4,"",6,"",8,"",10,15,20,25,30])

        ax_1.set_xticks(np.arange(0,V.shape[0]-20+1,20))
        # xticklabels = []
        # for yy in range(0,V.shape[0]-20+1,20):
        #     coord_str = str(valley_lat[yy]) + "N \n" + str(valley_lon[yy]) + "E"
        #     xticklabels.append(coord_str)

#               ax_1.set_xticklabels(xticklabels,rotation=-20,fontsize=10)

        for axi in [ax_1]:
            axi.fill_between(range(0,V.shape[0]),z[:,0,t],color="burlywood")
            #axi.plot(range(ddd,bl_height.shape[0]-ddd),bl_height[ddd:-ddd,t],":k",alpha=0.7)
            axi.set_facecolor("white")
            axi.grid(linestyle=":",alpha=0.3)
            axi.set_ylim(0,5000)
            axi.set_xlim(y_ll,y_ul)
            #ylabel = val_names[vv] + "\nHeight [m]"
            ylabel = "Height [m]"
            axi.set_title(val_name,fontsize=16)
            axi.set_ylabel(ylabel,fontsize=16)
            c=0
            for jj in range(0,4):
                axi.scatter(segs[vv,jj],profiles[0][int(segs[vv,jj])]-50,marker="x",s=140,color=colors_crosses[c])
                c+=1

            # axi.plot(ridges[0],color="k",linestyle=":")
            # axi.plot(ridges[1],color="k",linestyle=":")
            axi.plot(ridge_avg,color="k",linestyle="--",alpha=0.7)

        if vv==3: ax_1.set_xlabel("Along-valley grid points",fontsize=16)

    fname = "cross_sections_onepanel" + str(t) + ".pdf"
#    fname = "test_" + str(t) + ".pdf"
    plt.savefig(fname,bbox_inches = 'tight',dpi=300)
    #plt.savefig("test.pdf",bbox_inches = 'tight',dpi=300)
