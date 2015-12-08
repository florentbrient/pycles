import netCDF4 as nc
import numpy as np
import pylab as plt
import glob
import os





my_i = np.arange(700,720,1) # which profiles to include

n_i = np.shape(my_i)[0]
nkr = 32 #75
vars = ['temperature_mean','u_mean','v_mean','s_mean','viscosity_mean','diffusivity_mean' ,'wind_speed', 'wind_angle', 'theta_mean']
colors = ['White','White','LightPink','PowderBlue','HotPink','SkyBlue','MediumVioletRed','RoyalBlue','Maroon','Navy','White','Black']


dirs = glob.glob('./Output.Bomex.f827b*')


#--------MEANS
for dir in dirs:
    file = glob.glob(dir+'/stats/*')[0]
    data = nc.Dataset(file,'r')
    for var in vars:
        string = file[-5:-3]
        if string[0]=='_':
            string = string[1]
        print(string)
        if string == '3':
            pass
        else:
            title = var
            plt.figure(title)
            reg1 = np.zeros(nkr)
            for t in my_i:
                for k in np.arange(nkr):
                    reg1[k] += data.groups['profiles'].variables[var][t,k]/float(n_i)

            plt.plot(reg1, data.groups['profiles'].variables['z'][:], linewidth = 2,label = string, color=colors[int(string)])



for var in vars:
    plt.figure(var)
    plt.xlabel(var)
    plt.ylabel('z, m',fontsize=14 )
    plt.grid(True)
    plt.legend()


#----VARIANCES
variances = ['horizontal_velocity_variance','vertical_velocity_variance','resolved_tke','resolved_theta_variance']
for dir in dirs:
    file = glob.glob(dir+'/stats/*')[0]
    data = nc.Dataset(file,'r')
    string = file[-5:-3]
    if string[0]=='_':
        string = string[1]
    print(string)
    if string == '3' :
        pass
    else:
        reg1 = np.zeros(nkr)
        reg1w = np.zeros(nkr)
        regth = np.zeros(nkr)
        for t in my_i:
            for k in np.arange(nkr):
                reg1[k] += data.groups['profiles'].variables['u_mean2'][t,k]/float(n_i)-data.groups['profiles'].variables['u_mean'][t,k]**2/float(n_i)
                reg1[k] += data.groups['profiles'].variables['v_mean2'][t,k]/float(n_i)-data.groups['profiles'].variables['v_mean'][t,k]**2/float(n_i)
                reg1w[k] += data.groups['profiles'].variables['w_mean2'][t,k]/float(n_i)-data.groups['profiles'].variables['w_mean'][t,k]**2/float(n_i)
                regth[k] += data.groups['profiles'].variables['theta_mean2'][t,k]/float(n_i)-data.groups['profiles'].variables['theta_mean'][t,k]**2/float(n_i)

        plt.figure('horizontal_velocity_variance')
        plt.plot(reg1, data.groups['profiles'].variables['z'][:], linewidth = 2,label = string, color=colors[int(string)])

        plt.figure('vertical_velocity_variance')
        plt.plot(reg1w, data.groups['profiles'].variables['z'][:], linewidth = 2,label = string, color=colors[int(string)])


        plt.figure('resolved_tke')
        plt.plot((reg1[:]+reg1w[:])*0.5, data.groups['profiles'].variables['z'][:], linewidth = 2,label = string, color=colors[int(string)])

        plt.figure('resolved_theta_variance')
        plt.plot(regth, data.groups['profiles'].variables['z'][:], linewidth = 2,label = string, color=colors[int(string)])

for var in variances:
    plt.figure(var)
    plt.title(var)
    plt.ylabel('z, m',fontsize=14 )
    plt.xlabel(var, fontsize=14)
    plt.grid(True)
    plt.legend()

vars = ['b_flux_surface_mean', 'viscosity_max', 'viscosity_min','w_max', 'w_min', 'obukhov_length_mean','friction_velocity_mean']
for dir in dirs:
    file = glob.glob(dir+'/stats/*')[0]
    data = nc.Dataset(file,'r')
    for var in vars:
        string = file[-5:-3]
        if string[0]=='_':
            string = string[1]
        print(string)
        if string == '3':
            pass
        else:
            title = var
            plt.figure(title)
            plt.plot(data.groups['timeseries'].variables['t'][:],data.groups['timeseries'].variables[var][:],linewidth = 2,label = string, color=colors[int(string)])

for var in vars:
    plt.figure(var)
    plt.title(var)
    plt.ylabel(var,fontsize=14 )
    plt.xlabel('time', fontsize=14)
    plt.grid(True)
    plt.legend()

# vars = ['friction_velocity_mean']
#
#
# for var in vars:
#     nfig +=1
#     plt.figure(nfig)
#     plt.title(var+'_square')
#     plt.plot(smag.groups['timeseries'].variables['t'][:],smag.groups['timeseries'].variables[var][:]**2,'-b')
#     plt.plot(tke.groups['timeseries'].variables['t'][:],tke.groups['timeseries'].variables[var][:]**2,'-r')
#     plt.plot(m2.groups['timeseries'].variables['t'][:],m2.groups['timeseries'].variables[var][:]**2,'-g')
#
#
plt.show()
