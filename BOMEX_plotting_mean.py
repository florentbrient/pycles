import netCDF4 as nc
import numpy as np
import pylab as plt
import glob
import os





my_i = np.arange(180,217,1) # which profiles to include
#my_i = np.arange(1,2,1) # which profiles to include
print(my_i)
n_i = np.shape(my_i)[0]
nkr = 75
vars = ['cloud_fraction'] #,'u_mean','v_mean','s_mean','viscosity_mean','diffusivity_mean' ,'wind_speed', 'wind_angle', 'theta_mean']
colors = ['White','Black','LightPink','PowderBlue','HotPink','SkyBlue','MediumVioletRed','RoyalBlue','Maroon','Navy','White','Black']


dirs = glob.glob('./Output.Bomex.f827b*')


#--------MEANS
for dir in dirs:
    file = glob.glob(dir+'/stats/*')[0]
    data = nc.Dataset(file,'r')
    print(data.groups['profiles'].variables['cloud_fraction'][72,1:72])
    for var in vars:
        string = '1'
        if string[0]=='_':
            string = string[1]
#        print(string)
        if string == '3':
            pass
        else:
            title = var
            plt.figure(title)
            reg1 = np.zeros(nkr)
            for t in my_i:
            	#print(t)
                for k in np.arange(nkr):
                    #print(k)
                    reg1[k] += data.groups['profiles'].variables[var][t,k]/float(n_i)

            print(reg1)
            plt.plot(reg1, data.groups['profiles'].variables['z'][:], linewidth = 2,label = string, color='r')



for var in vars:
    plt.figure(var)
    #plt.xlabel(var)
    #plt.ylabel('z, m',fontsize=14 )
    plt.title('cloud_fraction_profiles',fontsize=14 )
    plt.grid(True)
    plt.legend()

plt.show()
