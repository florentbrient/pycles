import sys, os
sys.path.append('/Volumes/data5/ZTAN/PyRRTM/rrtmg_lw_v4.85/gcm_model_column_gfortran_f2py/build/f2py')

import read_input_lw
import rrtmg_lw
import numpy as np

read_input_lw.rrlw_con.pi = np.float64(np.arcsin(1)*2)

os.chdir('/Volumes/data5/ZTAN/PyRRTM/rrtmg_lw_v4.85/gcm_model_column_gfortran_f2py/runs_std_atm')

filename_in = './input_rrtm_MLS-cld-imca0' #x61
#filename_in = './input_rrtm_MLS-clr'     # x1
#filename_in = './input_rrtm_MLS-clr-aer12' # x51
#filename_in = './input_rrtm_MLS-clr-xsec' # x41
#filename_in = './input_rrtm_MLW-clr'     # x31
#filename_in = './input_rrtm_SAW-clr'     # x21
#filename_in = './input_rrtm_TROP-clr'    # x11

nlayers_out,iout_out,imca,icld_out,iaer_out,idrv,pavel_out,tavel_out,pz_out,tz_out,\
tbound_out,semiss_out,dtbound_out,coldry_out,wkl_out,wbrodl_out,wx_out,pwvcm_out,\
inflag_out,iceflag_out,liqflag_out,cldfrac_out,tauc,ciwp,clwp,rei,rel,tauaer_out = \
read_input_lw.read_input.readprof(filename=filename_in, 
                               cldfile='./IN_CLD_RRTM',
                               aerfile='./IN_AER_RRTM')
                               
rrtmg_lw.rrtmg_lw_init.rrtmg_lw_ini(np.float64(1.004e3))
ncol = 1
nlay = nlayers_out
icld = icld_out
play = np.copy(pavel_out[0:nlay], order='F').reshape((ncol,nlay))
plev = np.copy(pz_out[0:nlay+1], order='F').reshape((ncol,nlay+1))
tlay = np.copy(tavel_out[0:nlay], order='F').reshape((ncol,nlay))
tlev = np.copy(tz_out[0:nlay+1], order='F').reshape((ncol,nlay+1))
tsfc = np.copy(tbound_out, order='F').reshape(ncol)

h2ovmr = np.copy(wkl_out[0,0:nlay]/coldry_out[0:nlay], order='F').reshape((ncol,nlay))
co2vmr = np.copy(wkl_out[1,0:nlay]/coldry_out[0:nlay], order='F').reshape((ncol,nlay))
o3vmr  = np.copy(wkl_out[2,0:nlay]/coldry_out[0:nlay], order='F').reshape((ncol,nlay))
n2ovmr = np.copy(wkl_out[3,0:nlay]/coldry_out[0:nlay], order='F').reshape((ncol,nlay))
ch4vmr = np.copy(wkl_out[5,0:nlay]/coldry_out[0:nlay], order='F').reshape((ncol,nlay))
o2vmr  = np.copy(wkl_out[6,0:nlay]/coldry_out[0:nlay], order='F').reshape((ncol,nlay))

ccl4vmr = np.copy(wx_out[0,0:nlay]/coldry_out[0:nlay]*1.e+20, order='F').reshape((ncol,nlay))
cfc11vmr = np.copy(wx_out[1,0:nlay]/coldry_out[0:nlay]*1.e+20, order='F').reshape((ncol,nlay))
cfc12vmr = np.copy(wx_out[2,0:nlay]/coldry_out[0:nlay]*1.e+20, order='F').reshape((ncol,nlay))
cfc22vmr = np.copy(wx_out[3,0:nlay]/coldry_out[0:nlay]*1.e+20, order='F').reshape((ncol,nlay))
emis= np.copy(semiss_out, order='F').reshape((ncol,16))

inflglw  = inflag_out
iceflglw = iceflag_out
liqflglw = liqflag_out
cldfr = np.copy(cldfrac_out[0:nlay], order='F').reshape((ncol,nlay))

taucld = np.copy(tauc[:,0:nlay], order='F').reshape((16,ncol,nlay))
cicewp = np.copy(ciwp[0:nlay], order='F').reshape((ncol,nlay))
cliqwp = np.copy(clwp[0:nlay], order='F').reshape((ncol,nlay))
reice  = np.copy(rei[0:nlay], order='F').reshape((ncol,nlay))
reliq  = np.copy(rel[0:nlay], order='F').reshape((ncol,nlay))
tauaer = np.copy(tauaer_out[0:nlay,:], order='F').reshape((ncol,nlay,16))

uflx,dflx,hr,uflxc,dflxc,hrc,duflx_dt,duflxc_dt = rrtmg_lw.rrtmg_lw_rad.rrtmg_lw(
                                                     icld,idrv,play,plev,tlay,tlev,tsfc,
                                                     h2ovmr,o3vmr,co2vmr,ch4vmr,n2ovmr,o2vmr,
                                                     cfc11vmr,cfc12vmr,cfc22vmr,ccl4vmr,emis,
                                                     inflglw,iceflglw,liqflglw,cldfr,
                                                     taucld,cicewp,cliqwp,reice,reliq,tauaer,
                                                     ncol, nlay)
hr  = np.append(hr,0.0).reshape(ncol,nlay+1)
hrc = np.append(hrc,0.0).reshape(ncol,nlay+1)


print (' LEVEL    PRESSURE   UPWARD FLUX   DOWNWARD FLUX    NET FLUX       HEATING RATE')
print ('            mb          W/m2          W/m2           W/m2          degree/day')

for itemp in xrange(nlay+1):
    i = nlay-itemp
    if (pz_out[i]<1.0e-2):
        print (' %3d        %7.6f   %8.4f      %8.4f      %12.7f          %9.5f' % (i, pz_out[i],uflx[0,i],dflx[0,i],uflx[0,i]-dflx[0,i],hr[0,i]))
    elif (pz_out[i]<1.0e-1):
        print (' %3d         %6.5f    %8.4f      %8.4f      %12.7f          %9.5f' % (i, pz_out[i],uflx[0,i],dflx[0,i],uflx[0,i]-dflx[0,i],hr[0,i]))
    elif (pz_out[i]<1.0):
        print (' %3d        %6.4f     %8.4f      %8.4f      %12.7f          %9.5f' % (i, pz_out[i],uflx[0,i],dflx[0,i],uflx[0,i]-dflx[0,i],hr[0,i]))
    elif (pz_out[i]<10.0):
        print (' %3d       %6.3f      %8.4f      %8.4f      %12.7f          %9.5f' % (i, pz_out[i],uflx[0,i],dflx[0,i],uflx[0,i]-dflx[0,i],hr[0,i]))
    elif (pz_out[i]<100.0):
        print (' %3d      %6.2f       %8.4f      %8.4f      %12.7f          %9.5f' % (i, pz_out[i],uflx[0,i],dflx[0,i],uflx[0,i]-dflx[0,i],hr[0,i]))
    elif (pz_out[i]<1000.0):
        print (' %3d     %6.1f        %8.4f      %8.4f      %12.7f          %9.5f' % (i, pz_out[i],uflx[0,i],dflx[0,i],uflx[0,i]-dflx[0,i],hr[0,i]))
    else:
        print (' %3d     %6.1f        %8.4f      %8.4f      %12.7f          %9.5f' % (i, pz_out[i],uflx[0,i],dflx[0,i],uflx[0,i]-dflx[0,i],hr[0,i]))
