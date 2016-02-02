import sys, os
sys.path.append('/Volumes/data5/ZTAN/PyRRTM/rrtmg_sw_v3.9/gcm_model_column_gfortran_f2py/build/f2py')

import read_input_sw
import rrtmg_sw
import numpy as np

read_input_sw.rrsw_con.pi = np.float64(np.arcsin(1)*2)
rrsw_scon = 1.36822e+03

os.chdir('/Volumes/data5/ZTAN/PyRRTM/rrtmg_sw_v3.9/gcm_model_column_gfortran_f2py/runs_std_atm')

#filename_in = './input_rrtm_MLS-cld-imca0' #x61
filename_in = './input_rrtm_MLS-clr'     # x1
#filename_in = './input_rrtm_MLS-clr-aer12' # x51
#filename_in = './input_rrtm_MLS-clr-sza65' # x41
#filename_in = './input_rrtm_MLW-clr'     # x31
#filename_in = './input_rrtm_SAW-clr'     # x21
#filename_in = './input_rrtm_TROP-clr'    # x11

nlayers_out, iout_out, imca, icld_out, iaer_out, isccos_out, idelm_out, pdp, \
pavel_out, tavel_out, pz_out, tz_out, tbound_out, semiss_out, zenith_out, \
adjflux_out, dyofyr_out, adjes_out, coldry_out, wkl_out, inflag_out, iceflag_out,liqflag_out, \
cldfrac_out, tauc, ssac, asmc, fsfc, ciwp, clwp, rei, rel, tauaer_out, ssaaer_out, asmaer_out = \
read_input_sw.read_input.readprof(filename=filename_in, 
                               cldfile='./IN_CLD_RRTM',
                               aerfile='./IN_AER_RRTM')
                               
rrtmg_sw.rrtmg_sw_init.rrtmg_sw_ini(np.float64(1.004e3))
ncol = 1
nlay = nlayers_out
icld = icld_out
iaer = iaer_out

# p,t
play = np.copy(pavel_out[0:nlay], order='F').reshape((ncol,nlay))
plev = np.copy(pz_out[0:nlay+1], order='F').reshape((ncol,nlay+1))
tlay = np.copy(tavel_out[0:nlay], order='F').reshape((ncol,nlay))
tlev = np.copy(tz_out[0:nlay+1], order='F').reshape((ncol,nlay+1))
tsfc = np.copy(tbound_out, order='F').reshape(ncol)

# vmr
h2ovmr = np.copy(wkl_out[0,0:nlay]/coldry_out[0:nlay], order='F').reshape((ncol,nlay))
co2vmr = np.copy(wkl_out[1,0:nlay]/coldry_out[0:nlay], order='F').reshape((ncol,nlay))
o3vmr  = np.copy(wkl_out[2,0:nlay]/coldry_out[0:nlay], order='F').reshape((ncol,nlay))
n2ovmr = np.copy(wkl_out[3,0:nlay]/coldry_out[0:nlay], order='F').reshape((ncol,nlay))
ch4vmr = np.copy(wkl_out[5,0:nlay]/coldry_out[0:nlay], order='F').reshape((ncol,nlay))
o2vmr  = np.copy(wkl_out[6,0:nlay]/coldry_out[0:nlay], order='F').reshape((ncol,nlay))

# asdir, aldir, asdif, aldif <- semiss
asdir = 1.0 - np.copy(semiss_out[15], order='F').reshape(1) 
aldir = 1.0 - np.copy(semiss_out[15], order='F').reshape(1) 
asdif = 1.0 - np.copy(semiss_out[15], order='F').reshape(1) 
aldif = 1.0 - np.copy(semiss_out[15], order='F').reshape(1) 
# Note: albedo by band is not implemented !!

coszen = np.copy(zenith_out, order='F').reshape(1) 

inflgsw  = inflag_out
iceflgsw = iceflag_out
liqflgsw = liqflag_out
dyofyr   = dyofyr_out
adjes    = adjes_out
scon     = rrsw_scon 
idelm    = idelm_out

cldfr = np.copy(cldfrac_out[0:nlay], order='F').reshape((ncol,nlay))

taucld = np.copy(tauc[:,0:nlay], order='F').reshape((14,ncol,nlay))
ssacld = np.copy(ssac[:,0:nlay], order='F').reshape((14,ncol,nlay))
asmcld = np.copy(asmc[:,0:nlay], order='F').reshape((14,ncol,nlay))
fsfcld = np.copy(fsfc[:,0:nlay], order='F').reshape((14,ncol,nlay))

cicewp = np.copy(ciwp[0:nlay], order='F').reshape((ncol,nlay))
cliqwp = np.copy(clwp[0:nlay], order='F').reshape((ncol,nlay))
reice  = np.copy(rei[0:nlay], order='F').reshape((ncol,nlay))
reliq  = np.copy(rel[0:nlay], order='F').reshape((ncol,nlay))

tauaer = np.copy(tauaer_out[0:nlay,15:29], order='F').reshape((ncol,nlay,14))
ssaaer = np.copy(ssaaer_out[0:nlay,15:29], order='F').reshape((ncol,nlay,14))
asmaer = np.copy(asmaer_out[0:nlay,15:29], order='F').reshape((ncol,nlay,14))

ecaer  = np.zeros((ncol,nlay,6), dtype='float64', order='F') + 1.0e-15

# uflx,dflx,hr,uflxc,dflxc,hrc,duflx_dt,duflxc_dt = rrtmg_sw.rrtmg_sw_rad.rrtmg_sw(
#                                                     icld,idrv,play,plev,tlay,tlev,tsfc,
#                                                     h2ovmr,o3vmr,co2vmr,ch4vmr,n2ovmr,o2vmr,
#                                                     cfc11vmr,cfc12vmr,cfc22vmr,ccl4vmr,emis,
#                                                     inflgsw,iceflgsw,liqflgsw,cldfr,
#                                                     taucld,cicewp,cliqwp,reice,reliq,tauaer,
#                                                     ncol, nlay)

uflx, dflx, hr, uflxc, dflxc, hrc, dirdflx, difdflx = rrtmg_sw.rrtmg_sw_rad.rrtmg_sw(
                    icld, iaer, idelm, play,plev,tlay,tlev,tsfc,
                    h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr, 
                    asdir   ,asdif   ,aldir   ,aldif   , 
                    coszen  ,adjes   ,dyofyr  ,scon    ,
                    inflgsw, iceflgsw, liqflgsw, cldfr, 
                    taucld, ssacld, asmcld, fsfcld,
                    cicewp, cliqwp, reice, reliq,
                    tauaer, ssaaer, asmaer, ecaer,
                    ncol, nlay)  
                    
hr  = np.append(hr,0.0).reshape(ncol,nlay+1)
hrc = np.append(hrc,0.0).reshape(ncol,nlay+1)


print (' LEVEL PRESSURE   UPWARD FLUX   DIFDOWN FLUX  DIRDOWN FLUX  DOWNWARD FLUX   NET FLUX    HEATING RATE')
print ('         mb          W/m2          W/m2          W/m2        W/m2          W/m2       degree/day')

for itemp in xrange(nlay+1):
    i = nlay-itemp
    if (pz_out[i]<1.0e-2):
        print (' %3d   %7.6f   %10.4f    %10.4f    %10.4f    %10.4f    %11.6f    %10.5f' % (i, pz_out[i],uflx[0,i],difdflx[0,i],dirdflx[0,i],dflx[0,i],dflx[0,i]-uflx[0,i],hr[0,i]))
    elif (pz_out[i]<1.0e-1):
        print (' %3d    %6.5f   %10.4f    %10.4f    %10.4f    %10.4f    %11.6f    %10.5f' % (i, pz_out[i],uflx[0,i],difdflx[0,i],dirdflx[0,i],dflx[0,i],dflx[0,i]-uflx[0,i],hr[0,i]))
    elif (pz_out[i]<1.0):
        print (' %3d    %6.4f    %10.4f    %10.4f    %10.4f    %10.4f    %11.6f    %10.5f' % (i, pz_out[i],uflx[0,i],difdflx[0,i],dirdflx[0,i],dflx[0,i],dflx[0,i]-uflx[0,i],hr[0,i]))
    elif (pz_out[i]<10.0):
        print (' %3d    %6.3f    %10.4f    %10.4f    %10.4f    %10.4f    %11.6f    %10.5f' % (i, pz_out[i],uflx[0,i],difdflx[0,i],dirdflx[0,i],dflx[0,i],dflx[0,i]-uflx[0,i],hr[0,i]))
    elif (pz_out[i]<100.0):
        print (' %3d    %6.2f    %10.4f    %10.4f    %10.4f    %10.4f    %11.6f    %10.5f' % (i, pz_out[i],uflx[0,i],difdflx[0,i],dirdflx[0,i],dflx[0,i],dflx[0,i]-uflx[0,i],hr[0,i]))
    elif (pz_out[i]<1000.0):
        print (' %3d    %6.1f    %10.4f    %10.4f    %10.4f    %10.4f    %11.6f    %10.5f' % (i, pz_out[i],uflx[0,i],difdflx[0,i],dirdflx[0,i],dflx[0,i],dflx[0,i]-uflx[0,i],hr[0,i]))
    else:
        print (' %3d    %6.1f    %10.4f    %10.4f    %10.4f    %10.4f    %11.6f    %10.5f' % (i, pz_out[i],uflx[0,i],difdflx[0,i],dirdflx[0,i],dflx[0,i],dflx[0,i]-uflx[0,i],hr[0,i]))
