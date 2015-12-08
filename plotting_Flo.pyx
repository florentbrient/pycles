import netCDF4 as nc
import pylab as plt
import numpy as np
import matplotlib.cm as mplcm
import matplotlib.colors as colors

def main(files):

    nc_dataset_list(files)

    tmin = 300
    tmax = 649
    

    ps_vars = group_variables(files,'profiles')
    ts_vars = group_variables(files,'timeseries')

    for v in ps_vars:
        print v 
        try:
            plot_mean_profile(files,v, tmin, tmax)
        except:
            if v != 'z' and v != 't':
                plot_mean_profile(files,v, tmin, tmax)


    for v in ts_vars:
        try:
            plot_timeseries(files,v)
        except:
            pass



    close_datasets(files)
    return


def nc_dataset_list(files):
    files['datasets'] = []
    for p in files['paths']:
        files['datasets'].append(nc.Dataset(p,'r'))

    return

def close_datasets(files):
    for f in files['datasets']:
        f.close()

def return_data(files, group, varname):
    data = []
    for d in files['datasets']:
        data.append(d.groups[group].variables[varname])
    return data

def plot_timeseries(files,var_name):
    t = return_data(files, 'timeseries', 't')
    data = return_data(files, 'timeseries', var_name)

    NUM_COLORS = len(data)

    cm = plt.get_cmap('gist_rainbow')
    cNorm  = colors.Normalize(vmin=0, vmax=NUM_COLORS-1)
    scalarMap = mplcm.ScalarMappable(norm=cNorm, cmap=cm)
    
    fig = plt.figure(var_name + '_timeseries')
    ax = fig.add_subplot(111)
    ax.set_color_cycle([scalarMap.to_rgba(i) for i in range(NUM_COLORS)])
    count = 0
    for i in data:
        ax.plot(t[count][:]/3600.0,i[:],label=files['names'][count])
        count += 1
    ax.grid()
    ax.legend(loc=0)
    fig.savefig('./figs_new/' + var_name + '_timeseries' + '.png')
    plt.close()

    return

def group_variables(files, group):
    vars = files['datasets'][0].groups[group].variables.keys()
    return vars

def plot_mean_profile(files,var_name, tmin=0, tmax=-1):
    t = return_data(files, 'profiles', 't')
    z = return_data(files, 'profiles', 'z')
    data = return_data(files, 'profiles', var_name)
    
    NUM_COLORS = len(data)
    
    cm = plt.get_cmap('gist_rainbow')
    cNorm  = colors.Normalize(vmin=0, vmax=NUM_COLORS-1)
    scalarMap = mplcm.ScalarMappable(norm=cNorm, cmap=cm)
    
    fig = plt.figure(var_name + '_profiles')
    ax = fig.add_subplot(111)
    ax.set_color_cycle([scalarMap.to_rgba(i) for i in range(NUM_COLORS)])

    count = 0
    for i in data:
        ax.plot(np.mean(i[tmin:tmax,:],axis=0),z[count][:],label=files['names'][count])
        plt.title(var_name + '_profiles')
        count += 1
    ax.grid()
    ax.legend(loc=0)
    fig.savefig('./figs_new/' + var_name + '_profiles' + '.png')
    plt.close()

    return


def DYCOMS_files(files):

    files['paths'] = []
    files['names'] = []

    files['paths'].append('/Users/presselk/Desktop/DYCOMS/DYCOMS/Output.DYCOMS_RF01_order_11.d7b0a/stats/Stats.DYCOMS_RF01_order_11.nc')
    files['names'].append('WENO11')

    files['paths'].append('/Users/presselk/Desktop/DYCOMS/DYCOMS/Output.DYCOMS_RF01_order_9.31636/stats/Stats.DYCOMS_RF01_order_9.nc')
    files['names'].append('WENO9')

    files['paths'].append('/Users/presselk/Desktop/DYCOMS/DYCOMS/Output.DYCOMS_RF01_order_7.782e0/stats/Stats.DYCOMS_RF01_order_7.nc')
    files['names'].append('WENO7')

    files['paths'].append('/Users/presselk/Desktop/DYCOMS/DYCOMS/Output.DYCOMS_RF01_order_5.c7a5b/stats/Stats.DYCOMS_RF01_order_5.nc')
    files['names'].append('WENO5')

    files['paths'].append('/Users/presselk/Desktop/DYCOMS/DYCOMS_20m/Output.DYCOMS_RF01_order_3.a5558/stats/Stats.DYCOMS_RF01_order_3.nc')
    files['names'].append('WENO3')



    files['paths'].append('/Users/presselk/Desktop/DYCOMS/DYCOMS/Output.DYCOMS_RF01_order_2.ad106/stats/Stats.DYCOMS_RF01_order_2.nc')
    files['names'].append('Central 2nd')

    files['paths'].append('/Users/presselk/Desktop/DYCOMS/DYCOMS/Output.DYCOMS_RF01_order_4.87f0d/stats/Stats.DYCOMS_RF01_order_4.nc')
    files['names'].append('Central 4th')

    files['paths'].append('/Users/presselk/Desktop/DYCOMS/DYCOMS/Output.DYCOMS_RF01_order_6.f1591/stats/Stats.DYCOMS_RF01_order_6.nc')
    files['names'].append('Central 6th')

    files['paths'].append('/Users/presselk/Desktop/DYCOMS/DYCOMS/Output.DYCOMS_RF01_order_8.25d69/stats/Stats.DYCOMS_RF01_order_8.nc')
    files['names'].append('Central 8th')

    return

def DYCOMS_20m_files(files):

    files['paths'] = []
    files['names'] = []

    #files['paths'].append('./DYCOMS_20m/Output.DYCOMS_RF01_order_2.e403e/stats/Stats.DYCOMS_RF01_order_2.nc')
    #files['names'].append('Central 2nd')
    
    #files['paths'].append('./DYCOMS_20m/Output.DYCOMS_RF01_TKE_order_2.ee19f/stats/Stats.DYCOMS_RF01_TKE_order_2.nc')
    #files['names'].append('Central 2nd TKE')
    
    #files['paths'].append('./DYCOMS_10m/Output.DYCOMS_RF01_order_2_10m.29356/stats/Stats.DYCOMS_RF01_order_2_10m.nc')
    #files['names'].append('Central 2nd')


    #files['paths'].append('./DYCOMS_20m_TKE/Output.DYCOMS_RF01_order_2_TKE.cc614/stats/Stats.DYCOMS_RF01_order_2_TKE.nc')
    #files['names'].append('Central 2nd TKE')

    #files['paths'].append('./DYCOMS_20m/Output.DYCOMS_RF01_order_3.573c4/stats/Stats.DYCOMS_RF01_order_3.nc')
    #files['names'].append('WENO3')
    
    #files['paths'].append('./DYCOMS_20m/Output.DYCOMS_RF01_order_3mp.8d542/stats/Stats.DYCOMS_RF01_order_3mp.nc')
    #files['names'].append('WENO3mp')
    
    #files['paths'].append('./DYCOMS_20m/Output.DYCOMS_RF01_TKE_order_3.17945/stats/Stats.DYCOMS_RF01_TKE_order_3.nc')
    #files['names'].append('WENO3 TKE')

    files['paths'].append('./DYCOMS_20m/Output.DYCOMS_RF01_order_5.73a92/stats/Stats.DYCOMS_RF01_order_5.nc')
    files['names'].append('WENO5')
    
    files['paths'].append('./DYCOMS_20m/Output.DYCOMS_RF01_order_5mp.21f98/stats/Stats.DYCOMS_RF01_order_5mp.nc')
    files['names'].append('WENO5mp')
    
    
    #files['paths'].append('./DYCOMS_20m/Output.DYCOMS_RF01_TKE_order_5.7ec75/stats/Stats.DYCOMS_RF01_TKE_order_5.nc')
    #files['names'].append('WENO5 TKE')


    files['paths'].append('./DYCOMS_20m/Output.DYCOMS_RF01_order_7.5a01f/stats/Stats.DYCOMS_RF01_order_7.nc')
    files['names'].append('WENO7')
    
    files['paths'].append('./DYCOMS_20m/Output.DYCOMS_RF01_order_7mp.89cad/stats/Stats.DYCOMS_RF01_order_7mp.nc')
    files['names'].append('WENO7mp')

    files['paths'].append('./DYCOMS_20m/Output.DYCOMS_RF01_order_9.99393/stats/Stats.DYCOMS_RF01_order_9.nc')
    files['names'].append('WENO9')

    files['paths'].append('./DYCOMS_20m/Output.DYCOMS_RF01_order_9mp.14fef/stats/Stats.DYCOMS_RF01_order_9mp.nc')
    files['names'].append('WENO9mp')
    
    files['paths'].append('./DYCOMS_20m/Output.DYCOMS_RF01_order_11.759fd/stats/Stats.DYCOMS_RF01_order_11.nc')
    files['names'].append('WENO11')
    
    files['paths'].append('./DYCOMS_20m/Output.DYCOMS_RF01_order_11mp.8ea37/stats/Stats.DYCOMS_RF01_order_11mp.nc')
    files['names'].append('WENO11mp')

    #files['paths'].append('./DYCOMS_10m/Output.DYCOMS_RF01_order_7.8d25f/stats/Stats.DYCOMS_RF01_order_7.nc')
    #files['names'].append('WENO7')
    
    #files['paths'].append('./DYCOMS_20m/Output.DYCOMS_RF01_order_9.9b03a/stats/Stats.DYCOMS_RF01_order_9.nc')
    #files['names'].append('WENO9')
    
    #files['paths'].append('./DYCOMS_20m_TKE/Output.DYCOMS_RF01_order_9_TKE.336f8/stats/Stats.DYCOMS_RF01_order_9_TKE.nc')
    #files['names'].append('WENO9 TKE')
    
    #files['paths'].append('./DYCOMS_20m_TKE/Output.DYCOMS_RF01_order_9_TKE_pr.2458a/stats/Stats.DYCOMS_RF01_order_9_TKE_pr.nc')
    #files['names'].append('WENO9 TKE PR')
    
    #files['paths'].append('./DYCOMS_20m_TKE/Output.DYCOMS_RF01_order_11_TKE_pr.cc1ff/stats/Stats.DYCOMS_RF01_order_11_TKE_pr.nc')
    #files['names'].append('WENO11 TKE PR')
    

def DYCOMS_10m_files(files):

    files['paths'] = []
    files['names'] = []

    files['paths'].append('./DYCOMS_10m/Output.DYCOMS_RF01_order_2.55f1f/stats/Stats.DYCOMS_RF01_order_2.nc')
    files['names'].append('Central 2nd')

    files['paths'].append('./DYCOMS_10m/Output.DYCOMS_RF01_order_3.ff58c/stats/Stats.DYCOMS_RF01_order_3.nc')
    files['names'].append('WENO3')
    
    files['paths'].append('./DYCOMS_10m/Output.DYCOMS_RF01_order_5.5a92c/stats/Stats.DYCOMS_RF01_order_5.nc')
    files['names'].append('WENO5')
    
    files['paths'].append('./DYCOMS_10m/Output.DYCOMS_RF01_order_7.8b54f/stats/Stats.DYCOMS_RF01_order_7.nc')
    files['names'].append('WENO7')
    
    files['paths'].append('./DYCOMS_10m/Output.DYCOMS_RF01_order_9.020fb/stats/Stats.DYCOMS_RF01_order_9.nc')
    files['names'].append('WENO9')
    
    files['paths'].append('./DYCOMS_10m/Output.DYCOMS_RF01_order_11.533c6/stats/Stats.DYCOMS_RF01_order_11.nc')
    files['names'].append('WENO11')
    
def BOMEX(files): 
    files['paths'] = []
    files['names'] = []

    #files['paths'].append('./BOMEX/Output.Bomex_order_2.bee62/stats/Stats.Bomex_order_2.nc')
    #files['names'].append('Central 2nd')

    files['paths'].append('./BOMEX/Output.Bomex_order_3.5dc1f/stats/Stats.Bomex_order_3.nc')
    files['names'].append('WENO3')
    
    files['paths'].append('./BOMEX/Output.Bomex_order_4.decb3/stats/Stats.Bomex_order_4.nc')
    files['names'].append('Central 4th')
    
    files['paths'].append('./BOMEX/Output.Bomex_order_5.3a8a3/stats/Stats.Bomex_order_5.nc')
    files['names'].append('WENO5')
    
    files['paths'].append('./BOMEX/Output.Bomex_order_6.27f2e/stats/Stats.Bomex_order_6.nc')
    files['names'].append('Central 6th')

    files['paths'].append('./BOMEX/Output.Bomex_order_7.207ca/stats/Stats.Bomex_order_7.nc')
    files['names'].append('WENO7')
   
    files['paths'].append('./BOMEX/Output.Bomex_order_8.fcae6/stats/Stats.Bomex_order_8.nc')
    files['names'].append('Central 8th')
    
    files['paths'].append('./BOMEX/Output.Bomex_order_9.90a1b/stats/Stats.Bomex_order_9.nc')
    files['names'].append('WENO9')

    files['paths'].append('./BOMEX/Output.Bomex_order_11.219a7/stats/Stats.Bomex_order_11.nc')
    files['names'].append('WENO11')
    
def BOMEX_files(files): 
    files['paths'] = []
    files['names'] = []
    
    files['paths'].append('./Output.Var.bdbcb/stats/Stats.Var.nc')
    files['names'].append('Var +2 Bomex')
    
    files['paths'].append('./Output.Var.ed8da/stats/Stats.Var.nc')
    files['names'].append('Var -2 Bomex')
    
    
def VARF_files(files): 
    files['paths'] = []
    files['names'] = []

    files['paths'].append('./Output.VARF.ebde8/stats/Stats.VARF.nc')
    files['names'].append('Bomex ')
    
    files['paths'].append('./Output.VARF.965d3/stats/Stats.VARF.nc')
    files['names'].append('Bomex SH+0.25*SH')
    
    files['paths'].append('./Output.VARF.b91d0/stats/Stats.VARF.nc')
    files['names'].append('Bomex SH-0.25*SH')

if __name__ == '__main__':

    files = {}
    #BOMEX(files)
    #DYCOMS_10m_files(files)
    VARF_files(files)


    main(files)
