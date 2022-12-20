import numpy as np
import sys, os
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as colors
from matplotlib.cm import get_cmap
import netCDF4 as nc4


stations = ["pyramid","namche","lukla"]

titles_w = ["(a) NCO-P (27.959N 86.813E) 5050 masl","(c) Namche (27.802N 86.715E) 3570 masl","(e) Lukla (27.696N 86.723E) 2660 masl"]
titles_T = ["(b) NCO-P (27.959N 86.813E) 5050 masl","(d) Namche (27.802N 86.715E) 3570 masl","(f) Lukla (27.696N 86.723E) 2660 masl"]


xticks = range(0,5*24*2,12)
xlabels = ["06","12","18","00",
            "06","12","18","00",
            "06","12","18","00",
            "06","12","18","00",
            "06","12","18","00",]
xlabels_date = ["06","12\n17/12/2014","18","00",
            "06","12\n18/12/2014","18","00",
            "06","12\n19/12/2014","18","00",
            "06","12\n20/12/2014","18","00",
            "06","12\n21/12/2014","18","00",]

MAE_T_day = []
MAE_T_night = []
MAE_ws_day = []
MAE_ws_night = []

bias_T_day = []
bias_T_night = []
bias_ws_day = []
bias_ws_night = []

RMSE_T_day = []
RMSE_T_night = []
RMSE_ws_day = []
RMSE_ws_night = []

day_timesteps_30min = []
day_timesteps_60min = []
night_timesteps_30min = []
night_timesteps_60min = []

no_of_nans_day = []
no_of_nans_night = []

for dd in range(0,5): #separate day and night timesteps (5 days, 30min interval)
    for timestep in range(dd*24*2,dd*24*2+12*2):
        day_timesteps_30min.append(timestep)
    for timestep in range(dd*24,dd*24+12):
        day_timesteps_60min.append(timestep)
    for timestep in range(dd*24*2+12*2,dd*24*2+2*12*2):
        night_timesteps_30min.append(timestep)
    for timestep in range(dd*24+12,dd*24+2*12):
        night_timesteps_60min.append(timestep)

height_diff = [45.0,210.0,68.0] #height difference between the station and the model grid point
dry_lapserate = 9.8e-3
fig, axes = plt.subplots(3,2,figsize=(2*10,3*5))


for stat in range(3): #go through 3 stations
    stati = stations[stat]
    ax1 = axes[stat,0]
    ax2 = axes[stat,1]


    obs_ws = []
    infile = stati + "_ws.txt"
    with open(infile) as data_file:
        for line in data_file:
            parts = line.split() # split line into parts
            obs_ws.append(float(parts[0]))

    obs_ws = np.array(obs_ws)

    obs_wd = []
    infile = stati + "_wd.txt"
    with open(infile) as data_file:
        for line in data_file:
            parts = line.split() # split line into parts
            obs_wd.append(float(parts[0]))

    obs_temp = []
    infile = stati + "_temp.txt"
    with open(infile) as data_file:
        for line in data_file:
            parts = line.split() # split line into parts
            obs_temp.append(float(parts[0])+273.15)

    obs_temp = np.array(obs_temp)

    #read wrf data in, data extracted along the valley center lines
    data = stati + "_variables.nc"
    ncfile = nc4.Dataset(data)
    T2 =  np.array(ncfile.variables["T2"])
    U10 = np.array(ncfile.variables["U10"])
    V10 = np.array(ncfile.variables["V10"])

#    print(np.nanmean(T2))

    T2 = T2-dry_lapserate*height_diff[stat]

#    print(np.nanmean(T2))

    wspd10_wrf = np.sqrt(U10[0,:]**2 + V10[0,:]**2)

    wrf_Unorm = U10[0,:]/wspd10_wrf
    wrf_Vnorm = V10[0,:]/wspd10_wrf

    obs_Unorm = -np.sin(np.deg2rad(obs_wd[:]))
    obs_Vnorm = -np.cos(np.deg2rad(obs_wd[:]))

    obs_ws_GT1 = obs_ws.copy()
    obs_ws_LT1 = obs_ws.copy()


    #filter out timesteps with wind speed lower than anemometer treshold
    wind_treshold = 0.21 #m/s, http://www.proviento.com.pe/MW8008.CombiSD.pdf

    obs_ws_GT1[obs_ws_GT1 < wind_treshold] = np.NaN
    obs_ws_LT1[obs_ws_LT1 > wind_treshold] = np.NaN

    obs_Unorm[np.isnan(obs_ws_GT1)] = np.NaN
    obs_Vnorm[np.isnan(obs_ws_GT1)] = np.NaN

    no_of_nan_day = 0
    no_of_nan_night = 0

    #Wind speed 10m to 5m using Log-wind-profile
    z_0 = 0.5 #roughness length [m]
    log_fraction = (np.log(5.0/z_0)) / (np.log(10.0/z_0))
    wspd5_wrf = wspd10_wrf[:] * log_fraction

    print(np.amax(wspd5_wrf))

    if stat==0: #30min or 60min measurement frequency

        MAE_temp_day = []
        MAE_temp_night = []
        MAE_wspd_day = []
        MAE_wspd_night = []

        bias_temp_day = []
        bias_temp_night = []
        bias_wspd_day = []
        bias_wspd_night = []

        RMSE_temp_day = []
        RMSE_temp_night = []
        RMSE_wspd_day = []
        RMSE_wspd_night = []

        for ii in range(12*2,len(day_timesteps_30min)):
            bias_temp_day.append(T2[0,day_timesteps_30min[ii]]-obs_temp[day_timesteps_30min[ii]])
            bias_wspd_day.append(wspd5_wrf[day_timesteps_30min[ii]]-obs_ws_GT1[day_timesteps_30min[ii]])
            MAE_temp_day.append(abs(obs_temp[day_timesteps_30min[ii]]-T2[0,day_timesteps_30min[ii]]))
            MAE_wspd_day.append(abs(obs_ws_GT1[day_timesteps_30min[ii]]-wspd5_wrf[day_timesteps_30min[ii]]))
            RMSE_temp_day.append((abs(obs_temp[day_timesteps_30min[ii]]-T2[0,day_timesteps_30min[ii]]))**2)
            RMSE_wspd_day.append((abs(obs_ws_GT1[day_timesteps_30min[ii]]-wspd5_wrf[day_timesteps_30min[ii]]))**2)
            if(np.isnan(obs_ws_GT1[day_timesteps_30min[ii]])): no_of_nan_day += 1

        for ii in range(len(day_timesteps_30min)):
            bias_temp_night.append(T2[0,night_timesteps_30min[ii]]-obs_temp[night_timesteps_30min[ii]])
            bias_wspd_night.append(wspd5_wrf[night_timesteps_30min[ii]]-obs_ws_GT1[night_timesteps_30min[ii]])
            MAE_temp_night.append(abs(obs_temp[night_timesteps_30min[ii]]-T2[0,night_timesteps_30min[ii]]))
            MAE_wspd_night.append(abs(obs_ws_GT1[night_timesteps_30min[ii]]-wspd5_wrf[night_timesteps_30min[ii]]))
            RMSE_temp_night.append((abs(obs_temp[night_timesteps_30min[ii]]-T2[0,night_timesteps_30min[ii]]))**2)
            RMSE_wspd_night.append((abs(obs_ws_GT1[night_timesteps_30min[ii]]-wspd5_wrf[night_timesteps_30min[ii]]))**2)

            if(np.isnan(obs_ws_GT1[night_timesteps_30min[ii]])): no_of_nan_night += 1

        MAE_T_day.append(np.nanmean(MAE_temp_day))
        MAE_T_night.append(np.nanmean(MAE_temp_night))
        MAE_ws_day.append(np.nanmean(MAE_wspd_day))
        MAE_ws_night.append(np.nanmean(MAE_wspd_night))

        bias_T_day.append(np.nanmean(bias_temp_day))
        bias_T_night.append(np.nanmean(bias_temp_night))
        bias_ws_day.append(np.nanmean(bias_wspd_day))
        bias_ws_night.append(np.nanmean(bias_wspd_night))


        RMSE_T_day.append(np.sqrt(np.nanmean(RMSE_temp_day)))
        RMSE_T_night.append(np.sqrt(np.nanmean(RMSE_temp_night)))
        RMSE_ws_day.append(np.sqrt(np.nanmean(RMSE_wspd_day)))
        RMSE_ws_night.append(np.sqrt(np.nanmean(RMSE_wspd_night)))

    else:
        MAE_temp_day = []
        MAE_temp_night = []
        MAE_wspd_day = []
        MAE_wspd_night = []

        bias_temp_day = []
        bias_temp_night = []
        bias_wspd_day = []
        bias_wspd_night = []

        RMSE_temp_day = []
        RMSE_temp_night = []
        RMSE_wspd_day = []
        RMSE_wspd_night = []

        for ii in range(12,len(day_timesteps_60min)):
            bias_temp_day.append(T2[0,day_timesteps_30min[ii*2]]-obs_temp[day_timesteps_60min[ii]])
            bias_wspd_day.append(wspd5_wrf[day_timesteps_30min[ii*2]]-obs_ws_GT1[day_timesteps_60min[ii]])
            MAE_temp_day.append(abs(obs_temp[day_timesteps_60min[ii]]-T2[0,day_timesteps_30min[ii*2]]))
            MAE_wspd_day.append(abs(obs_ws_GT1[day_timesteps_60min[ii]]-wspd5_wrf[day_timesteps_30min[ii*2]]))
            RMSE_temp_day.append((abs(obs_temp[day_timesteps_60min[ii]]-T2[0,day_timesteps_30min[ii*2]]))**2)
            RMSE_wspd_day.append((abs(obs_ws_GT1[day_timesteps_60min[ii]]-wspd5_wrf[day_timesteps_30min[ii*2]]))**2)

            if(np.isnan(obs_ws_GT1[day_timesteps_60min[ii]])): no_of_nan_day += 1

        for ii in range(len(day_timesteps_60min)):
            bias_temp_night.append(T2[0,night_timesteps_30min[ii*2]]-obs_temp[night_timesteps_60min[ii]])
            bias_wspd_night.append(wspd5_wrf[night_timesteps_30min[ii*2]]-obs_ws_GT1[night_timesteps_60min[ii]])
            MAE_temp_night.append(abs(obs_temp[night_timesteps_60min[ii]]-T2[0,night_timesteps_30min[ii*2]]))
            MAE_wspd_night.append(abs(obs_ws_GT1[night_timesteps_60min[ii]]-wspd5_wrf[night_timesteps_30min[ii*2]]))
            RMSE_temp_night.append((abs(obs_temp[night_timesteps_60min[ii]]-T2[0,night_timesteps_30min[ii*2]]))**2)
            RMSE_wspd_night.append((abs(obs_ws_GT1[night_timesteps_60min[ii]]-wspd5_wrf[night_timesteps_30min[ii*2]]))**2)

            if(np.isnan(obs_ws_GT1[night_timesteps_60min[ii]])): no_of_nan_night += 1

        bias_T_day.append(np.nanmean(bias_temp_day))
        bias_T_night.append(np.nanmean(bias_temp_night))
        bias_ws_day.append(np.nanmean(bias_wspd_day))
        bias_ws_night.append(np.nanmean(bias_wspd_night))

        MAE_T_day.append(np.nanmean(MAE_temp_day))
        MAE_T_night.append(np.nanmean(MAE_temp_night))
        MAE_ws_day.append(np.nanmean(MAE_wspd_day))
        MAE_ws_night.append(np.nanmean(MAE_wspd_night))

        RMSE_T_day.append(np.sqrt(np.nanmean(RMSE_temp_day)))
        RMSE_T_night.append(np.sqrt(np.nanmean(RMSE_temp_night)))
        RMSE_ws_day.append(np.sqrt(np.nanmean(RMSE_wspd_day)))
        RMSE_ws_night.append(np.sqrt(np.nanmean(RMSE_wspd_night)))

    no_of_nans_day.append(no_of_nan_day)
    no_of_nans_night.append(no_of_nan_night)



    for i in range(0,5):
        for axi in [ax1,ax2]:
            axi.fill_between([i*48+24,(i+1)*48],[30,30],[-20,-20],color="gray",alpha=0.2)

    #ax1.plot(range(0,240),wspd_wrf,label="Model")

    ax1.set_ylabel("5 meter wind speed [m s$^{-1}$]")

    if stat==0:
#        ax1.plot(range(0,240),wspd10_wrf,label="Model 10m",color="black",linestyle=":")
        ax1.plot(range(0,240),wspd5_wrf,label="Model",color="black")
        ax1.plot(range(0,240),obs_ws,label="Observation",color="red")

        ax1.quiver(range(0,240)[::2],-1.5,wrf_Unorm[::2],wrf_Vnorm[::2],scale=50,width=0.002,pivot="middle",color="black")
        ax1.quiver(range(0,240)[::2],-3,obs_Unorm[::2],obs_Vnorm[::2],scale=50,width=0.002,pivot="middle",color="red")

        ax2.plot(range(0,240),obs_temp-273.15,color="red")
        ax2.plot(range(0,240),T2[0,:]-273.15,color="black")

    else:
#        ax1.plot(range(0,240),wspd10_wrf,label="Model 10m",color="black",linestyle=":")
        ax1.plot(range(0,240),wspd5_wrf,label="Model",color="black")
        ax1.plot(range(0,240)[::2],obs_ws_GT1,label="Observation",color="red",linestyle="-")
#        ax1.plot(range(0,240)[::2],obs_ws_LT1,label="Observation",color="orange")

        ax1.quiver(range(0,240)[::2],-1,wrf_Unorm[::2],wrf_Vnorm[::2],scale=50,width=0.002,pivot="middle",color="black")
        ax1.quiver(range(0,240)[::2],-2,obs_Unorm,obs_Vnorm,scale=50,width=0.002,pivot="middle",color="red")

        ax2.plot(range(0,240)[::2],obs_temp-273.15,color="red",)
        ax2.plot(range(0,240),T2[0,:]-273.15,color="black")

    ax1.plot([0,240],[0,0],color="gray")


    if stat in [0]: ax1.set_yticks([0,5,10,15,20])
    if stat in [1]: ax1.set_yticks([0,3,6,9,12])
    if stat in [2]: ax1.set_yticks([0,2,4,6,8])


    # ax1.set_yticks([0,2.5,5,7.5,10])
    # ax1.set_yticklabels(["0"," ","5"," ","10"])
#    ax1.set_ylim(-5,25)

    ax2.set_ylabel("2 meter temperature [$^\circ$C]")
#    ax2.set_ylim(-17,0)

    for axi in [ax1,ax2]:
        axi.set_xlim(24,240)
        axi.grid(linestyle="--",alpha=0.7)
        axi.set_xticks(xticks)
    ax1.set_title(titles_w[stat])
    ax2.set_title(titles_T[stat])



    if stat in [0]: ax1.set_ylim(-4,20)
    if stat in [1]: ax1.set_ylim(-2.5,12)
    if stat in [2]: ax1.set_ylim(-2.5,8)
    ax1.set_xticklabels([])
    ax2.set_xticklabels([])

print()
for stati in range(3):
    print(stations[stati])
    print("Bias T Daytime",np.round(bias_T_day[stati],1)," Kelvin")
    print("Bias T Night-time",np.round(bias_T_night[stati],1)," Kelvin")
    print("MAE T Daytime",np.round(MAE_T_day[stati],1)," Kelvin")
    print("MAE T Night-time",np.round(MAE_T_night[stati],1)," Kelvin")
    print("RMSE T Daytime",np.round(RMSE_T_day[stati],1)," Kelvin")
    print("RMSE T Night-time",np.round(RMSE_T_night[stati],1)," Kelvin")
    print("Bias WSPD Daytime",np.round(bias_ws_day[stati],1)," m s$^{-1}$","No. of NaN ",no_of_nans_day[stati])
    print("Bias WSPD Nighttime",np.round(bias_ws_night[stati],1)," m s$^{-1}$","No. of NaN ",no_of_nans_night[stati])
    print("MAE WSPD Daytime",np.round(MAE_ws_day[stati],1)," m s$^{-1}$","No. of NaN ",no_of_nans_day[stati])
    print("MAE WSPD Nighttime",np.round(MAE_ws_night[stati],1)," m s$^{-1}$","No. of NaN ",no_of_nans_night[stati])
    print("RMSE WSPD Daytime",np.round(RMSE_ws_day[stati],1)," m s$^{-1}$","No. of NaN ",no_of_nans_day[stati])
    print("RMSE WSPD Nighttime",np.round(RMSE_ws_night[stati],1)," m s$^{-1}$","No. of NaN ",no_of_nans_night[stati])
    print()



lgnd = axes[0,0].legend(bbox_to_anchor=(1.4, -2.55),ncol=2,fontsize=16)

#axes[0,0].text(-27.0,2.0,"5 meter wind [m s$^{-1}$]",color="red",rotation=90)
#axes[1,0].text(-23.5,0.75,"5 meter wind [m s$^{-1}$]",color="red",rotation=90)
#axes[2,0].text(-23.5,0.75,"5 meter wind [m s$^{-1}$]",color="red",rotation=90)

axes[0,1].set_ylim(-20,0)
axes[1,1].set_ylim(-10,10)
axes[2,1].set_ylim(-5,15)

axes[2,0].set_xticklabels(xlabels_date)
axes[2,1].set_xticklabels(xlabels_date)
axes[2,0].set_xlabel("Local time")
axes[2,1].set_xlabel("Local time")

for axis in axes:
    for axi in axis:
        axi.set_xlim(24,240)

#plt.savefig("test.pdf",bbox_inches = 'tight')
plt.savefig("obs_comparison_log.pdf",bbox_inches = 'tight',dpi=300)
