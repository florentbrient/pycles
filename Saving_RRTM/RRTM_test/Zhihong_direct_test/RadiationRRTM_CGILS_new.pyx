__author__ = 'pressel'
import numpy as np
cimport numpy as np
import sys
from libc.math cimport exp, fmin, fmax
from libc.stdlib cimport malloc, free
# from Parallel import get_z_pencils_double, z_pencils_inverse_double
# import pylab as plt
cimport ThermodynamicFunctions
import cython
from cython.parallel import prange
import ProfileStatistics
import timeit
from Parallel import reduce_to_root_2d

include 'Parameters.pxi'

from mpi4py import MPI

# Note: the RRTM modules are compiled in the 'RRTMG' directory:
cdef extern:
    void c_rrtmg_lw_init(double *cpdair)
    void c_rrtmg_lw (
             int *ncol    ,int *nlay    ,int *icld    ,int *idrv    , 
             double *play    ,double *plev    ,double *tlay    ,double *tlev    ,double *tsfc    , 
             double *h2ovmr  ,double *o3vmr   ,double *co2vmr  ,double *ch4vmr  ,double *n2ovmr  ,double *o2vmr, 
             double *cfc11vmr,double *cfc12vmr,double *cfc22vmr,double *ccl4vmr ,double *emis    , 
             int *inflglw ,int *iceflglw,int *liqflglw,double *cldfr   , 
             double *taucld  ,double *cicewp  ,double *cliqwp  ,double *reice   ,double *reliq   , 
             double *tauaer  , 
             double *uflx    ,double *dflx    ,double *hr      ,double *uflxc   ,double *dflxc,  double *hrc, 
             double *duflx_dt,double *duflxc_dt )
    void c_rrtmg_sw_init(double *cpdair)
    void c_rrtmg_sw (int *ncol    ,int *nlay    ,int *icld    ,int *iaer    , 
             double *play    ,double *plev    ,double *tlay    ,double *tlev    ,double *tsfc    , 
             double *h2ovmr  ,double *o3vmr   ,double *co2vmr  ,double *ch4vmr  ,double *n2ovmr  ,double *o2vmr, 
             double *asdir   ,double *asdif   ,double *aldir   ,double *aldif   , 
             double *coszen  ,double *adjes   ,int *dyofyr  ,double *scon    , 
             int *inflgsw ,int *iceflgsw,int *liqflgsw,double *cldfr   , 
             double *taucld  ,double *ssacld  ,double *asmcld  ,double *fsfcld  , 
             double *cicewp  ,double *cliqwp  ,double *reice   ,double *reliq   , 
             double *tauaer  ,double *ssaaer  ,double *asmaer  ,double *ecaer   , 
             double *swuflx  ,double *swdflx  ,double *swhr    ,double *swuflxc ,double *swdflxc ,double *swhrc)


import netCDF4 as nc

class RadiationRRTM_CGILS_new: # No change currently, should use the reference profile to extend the calculation.
    def __init__(self,grid):

        self.lw_flux_down = np.zeros(grid.nz)
        self.lw_flux_up = np.zeros(grid.nz)
        self.sw_flux_down = np.zeros(grid.nz)
        self.sw_flux_up = np.zeros(grid.nz)
        
        self.lw_flux_down_clear = np.zeros(grid.nz)
        self.lw_flux_up_clear = np.zeros(grid.nz)
        self.sw_flux_down_clear = np.zeros(grid.nz)
        self.sw_flux_up_clear = np.zeros(grid.nz)
        self.dTdt_clear = np.zeros(grid.nz)
        
        self.dTdt_rad = np.zeros((grid.nxl,grid.nyl,grid.nzl))
        self.dsdt_rad = np.zeros((grid.nxl,grid.nyl,grid.nzl))
        
        self.osr = np.zeros((grid.nxl, grid.nyl))
        self.olr = np.zeros((grid.nxl, grid.nyl))
        
        self.rel_tropo_T = None
        self.rel_toa_T   = None
        self.n_ext = None
        
        # CGILS-specific: 
        self.ref_aloft = None
        self.pressure_aloft = None
        self.p0i_aloft = None
        self.temperature_aloft = None
        self.y_vapor_aloft = None
        self.y_cond_aloft = None
        self.cloud_fraction_aloft = None
        
        # Test: 
        self.p0i_ext = None
        self.pressure_ext = None
        self.temperature_ext = None
        self.y_vapor_ext = None
        self.y_cond_ext = None
        self.cloud_fraction_ext = None
        
        self.is_vapor  = None
        self.is_liquid = None
        self.is_ice    = None
        
        # Added for PyRRTM
        self.lw_input_file = None
        self.lw_gas    = None
        self.trace     = None
        self.co2_factor = None
        self.h2o_factor = None
        self.lw_np = None
        self.uniform_reliq = None
        self.smooth_qc = None
        self.radiation_frequency = None
        
        self.o3vmr    = None
        self.co2vmr   = None
        self.ch4vmr   = None
        self.n2ovmr   = None
        self.o2vmr    = None
        self.cfc11vmr = None
        self.cfc12vmr = None
        self.cfc22vmr = None
        self.ccl4vmr  = None
        
        self.dyofyr   = None
        self.adjes    = None
        self.scon     = None
        self.coszen   = None
        self.adif     = None
        
        # Test
        self.play_in_ext  = None
        self.plev_in_ext  = None    
        self.tlay_in_ext  = None    
        self.tlev_in_ext  = None    
        self.h2ovmr_in_ext  = None
        self.cldfr_in_ext   = None
        self.cliqwp_in_ext  = None
        self.reliq_in_ext  = None
        self.tsfc_in_ext   = None
         
        self.debug_mode = None
        self.n_ext = None
        
        return

    def initialize(self,case_dict,grid,basicstate,scalars,thermodynamics,io):
    
        self.n_ext = 50
    
        try:
            self.ref_aloft = case_dict['radiation']['ref_aloft']
        except:
            print('Reference profile of radiation is not set, using default: ref_aloft = False')
            self.ref_aloft = False
            
        self.pressure_aloft    = np.zeros(self.n_ext,dtype=np.double,order='F')
        self.p0i_aloft         = np.zeros(self.n_ext+1 ,dtype=np.double,order='F')
        self.temperature_aloft = np.zeros(self.n_ext,dtype=np.double,order='F')
        self.y_vapor_aloft     = np.zeros(self.n_ext ,dtype=np.double,order='F')
        self.y_cond_aloft      = np.zeros(self.n_ext ,dtype=np.double,order='F')
        self.cloud_fraction_aloft = np.zeros(self.n_ext ,dtype=np.double,order='F')
            
        if self.ref_aloft:
            self.p0i_aloft[0] = basicstate.p0_i_global_ghosted[grid.gw + grid.nz-1]
            for k in xrange(1,self.n_ext+1):
                self.p0i_aloft[k] = self.p0i_aloft[k-1] - self.p0i_aloft[0]/self.n_ext
                self.pressure_aloft[k-1] = (self.p0i_aloft[k] + self.p0i_aloft[k-1])*0.5

            p = case_dict['forcing']['pressure']
            # Interpolate for temperature and moisture aloft (above of domain top) 
            try:
                self.temperature_aloft = np.interp(self.pressure_aloft,p,case_dict['forcing']['t_ref'])
            except:
                self.temperature_aloft = np.interp(self.pressure_aloft,p,case_dict['forcing']['thl_ref']*(p/basicstate.p00)**(Ra/cpa))
            self.y_vapor_aloft = np.interp(self.pressure_aloft,p,case_dict['forcing']['yv_ref'])    
            
        else: # ADDED BY ZTAN FOR EXTENDED TEMPERATURE PROFILE BEYOND LES TOP, ONLY USED WHEN self.ref_aloft = False
            try:
                self.rel_tropo_T = case_dict['forcing']['rel_tropo_T']
            except:
                print('Tropopause Temperature not set so RadiationFMS takes default value: rel_tropo_T = 240.0K .')
                self.rel_tropo_T = 240.0

            try:
                self.rel_toa_T = case_dict['forcing']['rel_toa_T']
            except:
                print('TOA Temperature not set so RadiationFMS takes default value: rel_toa_T = 220.0K .')
                self.rel_toa_T = 220.0
            
            try:
                self.rel_rh = case_dict['forcing']['rel_rh']
            except:
                print('Reference Relative Humidity not set so RadiationFMS takes default value: rel_rh = 0.4 .')
                self.rel_rh = 0.4
            
        # Added: calculate the extended pressure profiles
        pressure    = basicstate.p0_global_noghost
        p0i  = basicstate.p0_i_global_ghosted
        self.pressure_ext  = np.zeros(grid.nz + self.n_ext ,dtype=np.double,order='F')
        self.p0i_ext  = np.zeros(grid.nz + self.n_ext+1 ,dtype=np.double,order='F')
        
        for k in xrange(grid.nz + self.n_ext):
            if k <= grid.nz :
                self.p0i_ext[k] = p0i[grid.gw + k-1]
            elif self.ref_aloft:
                self.p0i_ext[k] = self.p0i_aloft[k-grid.nz]
            else:
                self.p0i_ext[k] = self.p0i_ext[k-1] - self.p0i_ext[grid.nz]/self.n_ext
        for k in xrange(grid.nz + self.n_ext):
            if k <= (grid.nz-1):
                self.pressure_ext[k] = pressure[k]
            else:
                self.pressure_ext[k] = (self.p0i_ext[k] + self.p0i_ext[k+1])*0.5 
                
            

        # Added: Initialize rrtmg_lw and rrtmg_sw
        cdef double cpdair = np.float64(cpa)
        c_rrtmg_lw_init(&cpdair)
        c_rrtmg_sw_init(&cpdair)
        
        # Added: Read in gas data (with CO2 multiplier)
        # Compute top of atmosphere solar insolation
        try:
            lw_input_file = case_dict['radiation']['lw_input_file']
        except:
            print('RRTM GHG values are not given a member of forcing case_dict dictionary. Required by RadiationFMS')
            sys.exit()
            
        try:
            self.co2_factor = case_dict['radiation']['co2_factor']
        except:
            print('CO2_factor not set so RadiationFMS takes default value: co2_factor = 1.0 .')
            self.co2_factor = 1.0  
            
        try:
            self.h2o_factor = case_dict['radiation']['h2o_factor']
        except:
            print('h2o_factor not set so RadiationFMS takes default value: h2o_factor = 1.0 .')
            self.h2o_factor = 1.0   
            
        try:
            self.uniform_reliq = case_dict['radiation']['uniform_reliq']
        except:
            print('uniform_reliq not set so RadiationFMS takes default value: uniform_reliq = False.')
            self.uniform_reliq = False  
        
        try:
            self.uniform_dTrad = case_dict['radiation']['uniform_dTrad']
        except:
            print('uniform_dTrad not set so RadiationFMS takes default value: uniform_dTrad = True.')
            self.uniform_dTrad = True
            
        try:
            self.smooth_qc = case_dict['radiation']['smooth_qc']
        except:
            # print('smooth_qc not set so RadiationFMS takes default value: smooth_qc = False.')
            self.smooth_qc = False
            
        try: # WRITE HERE !!
            o3_trace = case_dict['radiation']['o3mmr'][::-1]*28.97/47.9982   # O3 VMR (from SRF to TOP)
            o3_pressure = case_dict['forcing']['pressure'][::-1]/100.0       # Pressure (from SRF to TOP) in hPa
            # can't do simple interpolation... Need to conserve column path !!!
            self.use_o3in = True
        except:
            print('O3 profile not set so RadiationFMS takes default RRTM profile.')
            self.use_o3in = False
                   
        lw_gas = nc.Dataset(lw_input_file,  "r")
        
        # Trace Gases
        lw_pressure = np.asarray(lw_gas.variables['Pressure'])
        lw_absorber = np.asarray(lw_gas.variables['AbsorberAmountMLS'])
        lw_absorber = np.where(lw_absorber>2.0, np.zeros_like(lw_absorber), lw_absorber)
        lw_ngas = lw_absorber.shape[1]
        self.lw_np = lw_absorber.shape[0]

        # 9 Gases: O3, CO2, CH4, N2O, O2, CFC11, CFC12, CFC22, CCL4

        # From rad_driver.f90, lines 546 to 552
        trace = np.zeros((9,self.lw_np),dtype=np.double,order='F')
        for i in xrange(lw_ngas):
            gas_name = ''.join(lw_gas.variables['AbsorberNames'][i,:])
            if 'O3' in gas_name:
                trace[0,:] = lw_absorber[:,i].reshape(1,self.lw_np)
            elif 'CO2' in gas_name:
                trace[1,:] = lw_absorber[:,i].reshape(1,self.lw_np)*self.co2_factor
            elif 'CH4' in gas_name:
                trace[2,:] = lw_absorber[:,i].reshape(1,self.lw_np)
            elif 'N2O' in gas_name:
                trace[3,:] = lw_absorber[:,i].reshape(1,self.lw_np)
            elif 'O2' in gas_name:
                trace[4,:] = lw_absorber[:,i].reshape(1,self.lw_np)
            elif 'CFC11' in gas_name:
                trace[5,:] = lw_absorber[:,i].reshape(1,self.lw_np)
            elif 'CFC12' in gas_name:
                trace[6,:] = lw_absorber[:,i].reshape(1,self.lw_np)
            elif 'CFC22' in gas_name:
                trace[7,:] = lw_absorber[:,i].reshape(1,self.lw_np)
            elif 'CCL4' in gas_name:
                trace[8,:] = lw_absorber[:,i].reshape(1,self.lw_np)
                
        # From rad_driver.f90, lines 585 to 620       
        trpath = np.zeros((grid.nz + self.n_ext+1,9),dtype=np.double,order='F')
        plev = self.p0i_ext/100.0
        for i in xrange(1,grid.nz + self.n_ext+1):
            trpath[i,:] = trpath[i-1,:]
            if (plev[i-1] > lw_pressure[0]):
                trpath[i,:] = trpath[i,:] + (plev[i-1] - np.max((plev[i],lw_pressure[0])))/g*trace[:,0]
            for m in xrange(1,self.lw_np):
                #print i, m
                plow = np.min((plev[i-1],np.max((plev[i], lw_pressure[m-1]))))
                pupp = np.min((plev[i-1],np.max((plev[i], lw_pressure[m]))))
                if (plow > pupp):
                    pmid = 0.5*(plow+pupp)
                    wgtlow = (pmid-lw_pressure[m])/(lw_pressure[m-1]-lw_pressure[m])
                    wgtupp = (lw_pressure[m-1]-pmid)/(lw_pressure[m-1]-lw_pressure[m])
                    trpath[i,:] = trpath[i,:] + (plow-pupp)/g*(wgtlow*trace[:,m-1]  + wgtupp*trace[:,m])
            if (plev[i] < lw_pressure[self.lw_np-1]):
                trpath[i,:] = trpath[i,:] + (np.min((plev[i-1],lw_pressure[self.lw_np-1]))-plev[i])/g*trace[:,self.lw_np-1]

        tmpTrace = np.zeros((grid.nz + self.n_ext,9),dtype=np.double,order='F')
        for i in xrange(9):
            tmpTrace[:,i] = g/(plev[:-1]-plev[1:])*(trpath[1:,i]-trpath[:-1,i])

        self.o3vmr    = np.zeros(grid.nz + self.n_ext,dtype=np.double,order='F')
        self.co2vmr   = np.zeros(grid.nz + self.n_ext,dtype=np.double,order='F')
        self.ch4vmr   = np.zeros(grid.nz + self.n_ext,dtype=np.double,order='F')
        self.n2ovmr   = np.zeros(grid.nz + self.n_ext,dtype=np.double,order='F')
        self.o2vmr    = np.zeros(grid.nz + self.n_ext,dtype=np.double,order='F')
        self.cfc11vmr = np.zeros(grid.nz + self.n_ext,dtype=np.double,order='F')
        self.cfc12vmr = np.zeros(grid.nz + self.n_ext,dtype=np.double,order='F')
        self.cfc22vmr = np.zeros(grid.nz + self.n_ext,dtype=np.double,order='F')
        self.ccl4vmr  = np.zeros(grid.nz + self.n_ext,dtype=np.double,order='F')
        
        if self.use_o3in == False:
            self.o3vmr[:]  = tmpTrace[:,0]
        else:
            # o3_trace, o3_pressure
            trpath_o3 = np.zeros(grid.nz + self.n_ext+1,dtype=np.double,order='F')
            plev = self.p0i_ext/100.0
            self.o3_np = o3_trace.shape[0]
            for i in xrange(1,grid.nz + self.n_ext+1):
                trpath_o3[i] = trpath_o3[i-1]
                if (plev[i-1] > o3_pressure[0]):
                    trpath_o3[i] = trpath_o3[i] + (plev[i-1] - np.max((plev[i],o3_pressure[0])))/g*o3_trace[0]
                for m in xrange(1,self.o3_np):
                    #print i, m
                    plow = np.min((plev[i-1],np.max((plev[i], o3_pressure[m-1]))))
                    pupp = np.min((plev[i-1],np.max((plev[i], o3_pressure[m]))))
                    if (plow > pupp):
                        pmid = 0.5*(plow+pupp)
                        wgtlow = (pmid-o3_pressure[m])/(o3_pressure[m-1]-o3_pressure[m])
                        wgtupp = (o3_pressure[m-1]-pmid)/(o3_pressure[m-1]-o3_pressure[m])
                        trpath_o3[i] = trpath_o3[i] + (plow-pupp)/g*(wgtlow*o3_trace[m-1]  + wgtupp*o3_trace[m])
                if (plev[i] < o3_pressure[self.o3_np-1]):
                    trpath_o3[i] = trpath_o3[i] + (np.min((plev[i-1],o3_pressure[self.o3_np-1]))-plev[i])/g*o3_trace[self.o3_np-1]
            tmpTrace_o3 = np.zeros(grid.nz + self.n_ext,dtype=np.double,order='F')
            tmpTrace_o3[:] = g/(plev[:-1]-plev[1:])*(trpath_o3[1:]-trpath_o3[:-1])
            self.o3vmr[:] = tmpTrace_o3[:]
        
        self.co2vmr[:] = tmpTrace[:,1]
        self.ch4vmr[:] = tmpTrace[:,2]
        self.n2ovmr[:] = tmpTrace[:,3]
        self.o2vmr[:]  = tmpTrace[:,4]
        self.cfc11vmr[:] = tmpTrace[:,5]
        self.cfc12vmr[:] = tmpTrace[:,6]
        self.cfc22vmr[:] = tmpTrace[:,7]
        self.ccl4vmr[:]  = tmpTrace[:,8] 
          
        try:
            self.debug_mode = case_dict['radiation']['debug_mode']
        except:
            print('Debug Mode not set so RadiationFMS takes default value: debug_mode = False .')
            self.debug_mode = False
            
        try:
            self.dyofyr = case_dict['radiation']['dyofyr']
        except:
            print('Day of year not set so RadiationFMS takes default value: dyofyr = 0 .')
            self.dyofyr = 0
            
        try:
            self.adjes = case_dict['radiation']['adjes']
        except:
            print('Insolation adjustive factor not set so RadiationFMS takes default value: adjes = 0.5 (12 hour of daylight).')
            self.adjes = 0.5
            
        try:
            self.scon = case_dict['radiation']['solar_constant']
        except:
            print('Solar Constant not set so RadiationFMS takes default value: scon = 1360.0 .')
            self.scon = 1360.0
            
        try:
            self.coszen = case_dict['radiation']['coszen']
        except:
            print('Mean Daytime cos(SZA) not set so RadiationFMS takes default value: coszen = 2.0/pi .')
            self.coszen = 2.0/pi
            
        try:
            self.adif = case_dict['radiation']['adif']
        except:
            print('Surface diffusive albedo not set so RadiationFMS takes default value: adif = 0.06 .')
            self.adif = 0.06
            
        try:
            self.adir = case_dict['radiation']['adir']
        except:
            if (self.coszen > 0.0):
                self.adir = (.026/(self.coszen**1.7 + .065)+(.15*(self.coszen-0.10)*(self.coszen-0.50)*(self.coszen- 1.00)))
            else:
                self.adir = 0.0
            print('Surface direct albedo not set so RadiationFMS computes value: adif = %5.4f .'%(self.adir))
            
        try:
            self.radiation_frequency = case_dict['radiation']['radiation_frequency']
        except:
            print('radiation_frequency not set so RadiationFMS takes default value: radiation_frequency = 0.0 (compute at every step).')
            self.radiation_frequency = 0.0
            
        self.next_radiation_calculate = 0.0

        self.is_vapor  = True
        try:
            type(thermodynamics.y_vapor)
            print('RadiationFMS: Thermodynamics With Vapor')
        except:
            self.is_vapor = False
            print('RadiationFMS: Thermodynamics WithOUT Vapor')
            
        self.is_liquid  = True
        try:
            type(thermodynamics.y_liquid)
            print('RadiationFMS: Thermodynamics With Liquid')
        except:
            self.is_liquid = False   
            print('RadiationFMS: Thermodynamics WithOUT Liquid')
                     
        self.is_ice  = True
        try:
            type(thermodynamics.y_ice)
            print('RadiationFMS: Thermodynamics With Ice')
        except:
            self.is_ice = False
            print('RadiationFMS: Thermodynamics WithOUT Ice')
            
        #Initialize statistical output file
        self.init_output(grid, io)
            
        return

    
    @cython.boundscheck(False)
    @cython.wraparound(False)
    def update(self,grid,basicstate,scalars,velocities,thermodynamics,surface,comm,io,timestepping):
        # cdef double toc, tic
        cdef int gw = grid.gw
        cdef int nxl = grid.nxl
        cdef int nyl = grid.nyl
        cdef int nzl = grid.nzl
        
        cdef double [:,:,:] temperature = thermodynamics.temperature
        # Added liquid variables: will only work with (and only with) liquid !!
        cdef double [:,:,:] y_vapor  
        cdef double [:,:,:] y_liquid 
        cdef double [:,:,:] y_ice 
        
        cdef bint is_vapor  = self.is_vapor
        cdef bint is_liquid = self.is_liquid
        cdef bint is_ice    = self.is_ice
        
        if (self.is_vapor == True):
            y_vapor  = thermodynamics.y_vapor
            
        if (self.is_liquid == True):
            y_liquid = thermodynamics.y_liquid 
                       
        if (self.is_ice == True):
            y_ice    = thermodynamics.y_ice 
            
        cdef int i,j,k
        
        if (timestepping.rk_step==0): # Compute RRTM only on the first RK step
            if self.radiation_frequency <= 0.0:  # Frequency = 0.0: RRTM every step
                self.update_calculate(grid,basicstate,scalars,velocities,thermodynamics,surface,comm,io,timestepping)
            elif timestepping.time >= self.next_radiation_calculate: # Frequency > 0.0: RRTM whenever time > n*freq
                self.update_calculate(grid,basicstate,scalars,velocities,thermodynamics,surface,comm,io,timestepping)
                self.next_radiation_calculate = (timestepping.time//self.radiation_frequency + 1.0) * self.radiation_frequency
        
        # tic=timeit.default_timer()
        cdef double [:] total_tt_ext = self.total_tt_ext
        cdef double [:,:,:] total_tt_ext_3d = self.total_tt_ext_3d # Added for 3D output
        cdef bint uniform_dTrad = self.uniform_dTrad
        
        #Now add radiative temperature tendency
        cdef int y_water_dof = scalars.get_dof('y_water')
        cdef int specific_entropy_dof = scalars.get_dof('specific_entropy')
        cdef double [:,:,:] y_water = scalars.values[:,:,:,y_water_dof]
        cdef double [:,:,:] specific_entropy = scalars.values[:,:,:,specific_entropy_dof]
        cdef double [:,:,:] specific_entropy_t = scalars.tendencies[:,:,:,specific_entropy_dof]
        # cdef double [:,:,:] y_vapor = thermodynamics.y_vapor
        # cdef double [:,:,:] y_liquid = thermodynamics.y_liquid
        cdef double [:] alpha0 = basicstate.alpha0
        cdef double [:] p0 = basicstate.p0
        
        cdef double [:,:,:] tendency_T = self.dTdt_rad[:,:,:]
        cdef double [:,:,:] tendency_s = self.dsdt_rad[:,:,:]

        with nogil:
            for i in prange(gw,nxl-gw,schedule='static'):
                for j in xrange(gw,nyl-gw):
                    for k in xrange(gw,nzl-gw):
                        if uniform_dTrad:
                            tendency_T [i,j,k] = total_tt_ext[k-gw] #total_tt[k-gw]
                        else:
                            tendency_T [i,j,k] = total_tt_ext_3d[i-gw, j-gw, k-gw]
                        
                        if is_vapor:
                            if is_liquid:
                                if is_ice: # With Vapor, Liquid and Ice
                                    specific_entropy_t[i,j,k] = specific_entropy_t[i,j,k] + ( ThermodynamicFunctions.entropy_tendency(temperature[i,j,k],
                                                                                             p0[k],alpha0[k], y_vapor[i,j,k],y_liquid[i,j,k],
                                                                                             y_ice[i,j,k], 0.0,0.0,0.0,tendency_T [i,j,k]) )
                                    tendency_s [i,j,k] = ( ThermodynamicFunctions.entropy_tendency(temperature[i,j,k],
                                                          p0[k],alpha0[k], y_vapor[i,j,k],y_liquid[i,j,k],
                                                          y_ice[i,j,k], 0.0,0.0,0.0,tendency_T [i,j,k]) )
                                else: # With Vapor and Liquid, No Ice
                                    specific_entropy_t[i,j,k] = specific_entropy_t[i,j,k] + ( ThermodynamicFunctions.entropy_tendency(temperature[i,j,k],
                                                                                             p0[k],alpha0[k], y_vapor[i,j,k],y_liquid[i,j,k],
                                                                                             0.0, 0.0,0.0,0.0,tendency_T [i,j,k]) )
                                    tendency_s [i,j,k] = ( ThermodynamicFunctions.entropy_tendency(temperature[i,j,k],
                                                          p0[k],alpha0[k], y_vapor[i,j,k],y_liquid[i,j,k],
                                                          0.0, 0.0,0.0,0.0,tendency_T [i,j,k]) )
                            else: # With Vapor, No Liquid or Ice
                                specific_entropy_t[i,j,k] = specific_entropy_t[i,j,k] + ( ThermodynamicFunctions.entropy_tendency(temperature[i,j,k],
                                                                                         p0[k],alpha0[k], y_vapor[i,j,k],0.0,
                                                                                         0.0, 0.0,0.0,0.0,tendency_T [i,j,k]) )
                                tendency_s [i,j,k] = ( ThermodynamicFunctions.entropy_tendency(temperature[i,j,k],
                                                      p0[k],alpha0[k], y_vapor[i,j,k],0.0,
                                                      0.0, 0.0,0.0,0.0,tendency_T [i,j,k]) )
                        else: # No Vapor, Liquid or Ice
                            specific_entropy_t[i,j,k] = specific_entropy_t[i,j,k] + ( ThermodynamicFunctions.entropy_tendency(temperature[i,j,k],
                                                                                     p0[k],alpha0[k], 0.0,0.0,
                                                                                     0.0, 0.0,0.0,0.0,tendency_T [i,j,k]) )
                            tendency_s [i,j,k] = ( ThermodynamicFunctions.entropy_tendency(temperature[i,j,k],
                                                  p0[k],alpha0[k], 0.0,0.0,
                                                  0.0, 0.0,0.0,0.0,tendency_T [i,j,k]) )  
        #if comm.rank==0:
        #time.sleep(comm.rank/10.0)
        #print 'rk step is: ', timestepping.rk_step
        #print 'comm.rank is: ', comm.rank
        #print 'dTdt_rad, mean', np.mean(np.mean(self.dTdt_rad, axis=1),axis=0)
        #print 'dsdt_rad, mean', np.mean(np.mean(self.dsdt_rad, axis=1),axis=0)
        #print 'dsdt_all, mean', np.mean(np.mean(scalars.tendencies[:,:,:,specific_entropy_dof], axis=1),axis=0)
        #toc=timeit.default_timer()
        # print 'Timer 10: ', toc-tic
        #tic=toc
    
    

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def update_calculate(self,grid,basicstate,scalars,velocities,thermodynamics,surface,comm,io,timestepping):


        cdef double toc, tic
        # tic=timeit.default_timer() 
         
        cdef int gw = grid.gw
        cdef int nzg = grid.nzg
        cdef int nxl = grid.nxl
        cdef int nyl = grid.nyl
        cdef int nzl = grid.nzl
        cdef int kmin = grid.global_kmin
        cdef int nz = grid.nz
        cdef int nx = grid.nx
        cdef int ny = grid.ny
        cdef int ishift, jshift, kshift, ijshift
        cdef int n_ext = self.n_ext

        cdef double [:,:,:] temperature = thermodynamics.temperature
        # Added liquid variables: will only work with (and only with) liquid !!
        cdef double [:,:,:] y_vapor  
        cdef double [:,:,:] y_liquid 
        cdef double [:,:,:] y_ice 
        
        cdef bint is_vapor  = self.is_vapor
        cdef bint is_liquid = self.is_liquid
        cdef bint is_ice    = self.is_ice
        
        cdef bint uniform_reliq = self.uniform_reliq
        cdef bint smooth_qc = self.smooth_qc
        
        if (self.is_vapor == True):
            y_vapor  = thermodynamics.y_vapor
            
        if (self.is_liquid == True):
            y_liquid = thermodynamics.y_liquid 
                       
        if (self.is_ice == True):
            y_ice    = thermodynamics.y_ice 
        
        # LOOK FROM HERE ON ...
        # Don't average T, q, cc here - rather, do it later
        # Note: this code does not work with 3D-decomposition currently
        
        cdef int i,j,k
        cdef double [:] temperature_mean_ = np.zeros(nz,dtype=np.double,order='F')
        cdef double [:,:,:] temperature_ext = np.zeros((nxl-2*gw,nyl-2*gw,nz + n_ext),dtype=np.double,order='F')
        cdef double [:,:,:] y_vapor_ext     = np.zeros((nxl-2*gw,nyl-2*gw,nz + n_ext),dtype=np.double,order='F')
        cdef double [:,:,:] y_cond_ext      = np.zeros((nxl-2*gw,nyl-2*gw,nz + n_ext),dtype=np.double,order='F')
        cdef double [:,:,:] cloud_fraction_ext = np.zeros((nxl-2*gw,nyl-2*gw,nz + n_ext),dtype=np.double,order='F')
        
        # Added for 3D output
        cdef double [:,:,:] total_tt_ext_3d = np.zeros((nxl-2*gw,nyl-2*gw,nz + n_ext),dtype=np.double,order='F')
        
        # toc=timeit.default_timer()
        # print 'Timer 1: ', toc-tic
        # tic=toc
        
        with nogil:
            for i in xrange(nxl-2*gw):
                for j in xrange(nyl-2*gw):
                    for k in xrange(nzl-2*gw):
                        kshift = k + gw
                        ishift = i + gw
                        jshift = j + gw
                        temperature_mean_[kmin + k] += temperature[ishift,jshift,kshift]
                        temperature_ext[i,j,kmin + k] = temperature[ishift,jshift,kshift]
                        
                        # Need water vapor, liquid water, cloud fraction too...
                        if is_vapor:
                            y_vapor_ext[i,j,kmin + k] = y_vapor[ishift,jshift,kshift]
                            
                            if smooth_qc:  # Added 09/15/2015: smoothing
                                if is_liquid:
                                    if is_ice: # With Vapor, Liquid and Ice
                                        y_cond_ext[i,j,kmin + k] =  (     (y_liquid[ishift,jshift,kshift-1] + y_ice[ishift,jshift,kshift-1]) + 
                                                                      2.0*(y_liquid[ishift,jshift,kshift]   + y_ice[ishift,jshift,kshift])  + 
                                                                          (y_liquid[ishift,jshift,kshift+1] + y_ice[ishift,jshift,kshift+1]) ) / 4.0
                                        
                                    else: # With Vapor and Liquid, No Ice
                                        y_cond_ext[i,j,kmin + k] =  (     y_liquid[ishift,jshift,kshift-1] + 
                                                                      2.0*y_liquid[ishift,jshift,kshift]   + 
                                                                          y_liquid[ishift,jshift,kshift+1] ) / 4.0
                                else: # With Vapor, No Liquid or Ice
                                    y_cond_ext[i,j,kmin + k] =  (     y_ice[ishift,jshift,kshift-1] + 
                                                                  2.0*y_ice[ishift,jshift,kshift]   + 
                                                                      y_ice[ishift,jshift,kshift+1] ) / 4.0
                                    
                            else: 
                                if is_liquid:
                                    if is_ice: # With Vapor, Liquid and Ice
                                        y_cond_ext[i,j,kmin + k] = (y_liquid[ishift,jshift,kshift] + y_ice[ishift,jshift,kshift])
                                    else: # With Vapor and Liquid, No Ice
                                        y_cond_ext[i,j,kmin + k] =  y_liquid[ishift,jshift,kshift]
                                else: # With Vapor, No Liquid or Ice
                                    y_cond_ext[i,j,kmin + k] = y_ice[ishift,jshift,kshift]
                                
                                
                            if (y_cond_ext[i,j,kmin + k] > 1.0e-8):
                                cloud_fraction_ext[i,j,kmin + k] = 1.0  
                                                        
                        # else: # No Vapor, Liquid or Ice


        # toc=timeit.default_timer()
        # print 'Timer 2: ', toc-tic
        # tic=toc
        
        cdef double [:] temperature_mean = np.empty(nz,dtype=np.double,order='F')
        comm.cart_comm.Allreduce(temperature_mean_, temperature_mean,op = MPI.SUM)
        
        cdef double nxny = np.double(nx * ny)
        for k in xrange(nz):
            temperature_mean[k] = temperature_mean[k]/nxny


        # cdef double [:] temperature_ext = np.zeros(nz + n_ext ,dtype=np.double,order='F')
        cdef double [:] pressure      = basicstate.p0_global_noghost
        cdef double [:] pressure_ext  = self.pressure_ext
        
        cdef double [:] p0i_ext  = self.p0i_ext
        cdef double [:] p0i  = basicstate.p0_i_global_ghosted
        
        cdef bint ref_aloft = self.ref_aloft
        cdef double [:] temperature_aloft    = self.temperature_aloft
        cdef double [:] y_vapor_aloft        = self.y_vapor_aloft
        cdef double [:] y_cond_aloft         = self.y_cond_aloft
        cdef double [:] cloud_fraction_aloft = self.cloud_fraction_aloft

        cdef ThermodynamicFunctions.ref_adiabat_struct _ref
        cdef double t00   = temperature_mean[nz-1] 
        cdef double p00   = pressure[nz-1]  #pressure_ext[nz-1]
        
        # Compute the extended T-profile from a saturated adiabat 
        # (cut off at rel_tropo_T and linearly extend to rel_toa_T)
        # q-profile: Assume RH to be rel_rh
        
        if ref_aloft:
            with nogil:
                for k in xrange(nz, nz + n_ext):
                    for i in xrange(nxl-2*gw):
                        for j in xrange(nyl-2*gw):
                            temperature_ext[i,j,k]    = temperature_aloft[k-nz]
                            y_vapor_ext[i,j,k]        = y_vapor_aloft[k-nz]
                            y_cond_ext[i,j,k]         = y_cond_aloft[k-nz]   
                            cloud_fraction_ext[i,j,k] = cloud_fraction_aloft[k-nz] 
        else:
            _ref = ThermodynamicFunctions.compute_adiabat(p00, t00, pressure_ext[nz:], 1.0, self.rel_rh, 1.0, 
                    self.rel_tropo_T, self.rel_toa_T, n_ext) 
            with nogil:
                for k in xrange(nz, nz + n_ext):
                    for i in xrange(nxl-2*gw):
                        for j in xrange(nyl-2*gw):
                            temperature_ext[i,j,k] = _ref.t_ref[k-nz]
                            y_vapor_ext[i,j,k] = _ref.yv1_ref[k-nz]
                     
            free(_ref.t_ref )
            free(_ref.yv1_ref)
            free(_ref.yv2_ref)        
                
        cdef double h2o_factor = self.h2o_factor
        # toc=timeit.default_timer()
        # print 'Timer 3: ', toc-tic
        # tic=toc
                
        # Create input for rrtmg_lw and rrtmg_sw modules and call modules
        # LW GHG (size: nz+n_ext)
        cdef double surface_temperature = surface.sst
        cdef double [:] o3vmr   = self.o3vmr
        cdef double [:] co2vmr  = self.co2vmr
        cdef double [:] ch4vmr  = self.ch4vmr
        cdef double [:] n2ovmr  = self.n2ovmr
        cdef double [:] o2vmr   = self.o2vmr
        cdef double [:] cfc11vmr  = self.cfc11vmr
        cdef double [:] cfc12vmr  = self.cfc12vmr
        cdef double [:] cfc22vmr  = self.cfc22vmr
        cdef double [:] ccl4vmr   = self.ccl4vmr
        
        # Now we have 3d arrays: temperature_ext, y_vapor_ext, y_cond_ext, cloud_fraction_ext (nxl-2*gw, nyl-2*gw, nz+n_ext)
        # And 1d arrays: o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr, cfc11vmr, cfc12vmr, cfc22vmr, ccl4vmr (n_ext)
        # First: create input arrays for rrtmg: play, plev, tlay, tlev, tsfc, ... for rrtmg
        # should loop through all points
        # Second: call both subroutines
        # Third: convert the output arrays to the output array size
        
        cdef:
            int ncol = (nxl-2*gw)*(nyl-2*gw)
            int nlay = nz+n_ext
            int icld = 1
            int idrv = 0
            int iaer = 0
            int inflglw = 2
            int iceflglw = 3
            int liqflglw = 1
            int inflgsw = 2
            int iceflgsw = 3
            int liqflgsw = 1
            double [:,:] play_in = np.zeros((ncol,nz + n_ext),dtype=np.double,order='F')
            double [:,:] plev_in = np.zeros((ncol,nz + n_ext +1),dtype=np.double,order='F')
            double [:,:] tlay_in = np.zeros((ncol,nz + n_ext),dtype=np.double,order='F')
            double [:,:] tlev_in = np.zeros((ncol,nz + n_ext +1),dtype=np.double,order='F')
            double [:] tsfc_in = np.zeros((ncol),dtype=np.double,order='F')
            double [:,:] h2ovmr_in = np.zeros((ncol,nz + n_ext),dtype=np.double,order='F')
            double [:,:] o3vmr_in  = np.zeros((ncol,nz + n_ext),dtype=np.double,order='F')
            double [:,:] co2vmr_in = np.zeros((ncol,nz + n_ext),dtype=np.double,order='F')
            double [:,:] ch4vmr_in = np.zeros((ncol,nz + n_ext),dtype=np.double,order='F')
            double [:,:] n2ovmr_in = np.zeros((ncol,nz + n_ext),dtype=np.double,order='F')
            double [:,:] o2vmr_in  = np.zeros((ncol,nz + n_ext),dtype=np.double,order='F')
            double [:,:] cfc11vmr_in = np.zeros((ncol,nz + n_ext),dtype=np.double,order='F')
            double [:,:] cfc12vmr_in = np.zeros((ncol,nz + n_ext),dtype=np.double,order='F')
            double [:,:] cfc22vmr_in = np.zeros((ncol,nz + n_ext),dtype=np.double,order='F')
            double [:,:] ccl4vmr_in = np.zeros((ncol,nz + n_ext),dtype=np.double,order='F')
            double [:,:] emis_in = np.zeros((ncol,16),dtype=np.double,order='F')
            double [:,:] cldfr_in  = np.zeros((ncol,nz + n_ext),dtype=np.double,order='F')
            double [:,:] cicewp_in = np.zeros((ncol,nz + n_ext),dtype=np.double,order='F')
            double [:,:] cliqwp_in = np.zeros((ncol,nz + n_ext),dtype=np.double,order='F')
            double [:,:] reice_in  = np.zeros((ncol,nz + n_ext),dtype=np.double,order='F')
            double [:,:] reliq_in  = np.zeros((ncol,nz + n_ext),dtype=np.double,order='F')
            int dyofyr_in = self.dyofyr
            double  adjes_in  = self.adjes
            double  scon_in   = self.scon
            double  coszen    = self.coszen
            double  adif      = self.adif
            double  adir      = self.adir
            double [:] coszen_in = np.zeros((ncol),dtype=np.double,order='F')
            double [:] asdir_in = np.zeros((ncol),dtype=np.double,order='F')
            double [:] asdif_in = np.zeros((ncol),dtype=np.double,order='F')
            double [:] aldir_in = np.zeros((ncol),dtype=np.double,order='F')
            double [:] aldif_in = np.zeros((ncol),dtype=np.double,order='F')
            double [:,:,:] taucld_lw_in  = np.zeros((16,ncol,nz + n_ext),dtype=np.double,order='F')
            double [:,:,:] tauaer_lw_in  = np.zeros((ncol,nz + n_ext,16),dtype=np.double,order='F')
            double [:,:,:] taucld_sw_in  = np.zeros((14,ncol,nz + n_ext),dtype=np.double,order='F')
            double [:,:,:] ssacld_sw_in  = np.zeros((14,ncol,nz + n_ext),dtype=np.double,order='F')
            double [:,:,:] asmcld_sw_in  = np.zeros((14,ncol,nz + n_ext),dtype=np.double,order='F')
            double [:,:,:] fsfcld_sw_in  = np.zeros((14,ncol,nz + n_ext),dtype=np.double,order='F')
            double [:,:,:] tauaer_sw_in  = np.zeros((ncol,nz + n_ext,14),dtype=np.double,order='F')
            double [:,:,:] ssaaer_sw_in  = np.zeros((ncol,nz + n_ext,14),dtype=np.double,order='F')
            double [:,:,:] asmaer_sw_in  = np.zeros((ncol,nz + n_ext,14),dtype=np.double,order='F')
            double [:,:,:] ecaer_sw_in  = np.zeros((ncol,nz + n_ext,6),dtype=np.double,order='F')
            
            # Output
            double[:,:] uflx_lw_out = np.zeros((ncol,nz + n_ext +1),dtype=np.double,order='F')
            double[:,:] dflx_lw_out = np.zeros((ncol,nz + n_ext +1),dtype=np.double,order='F')
            double[:,:] hr_lw_out = np.zeros((ncol,nz + n_ext),dtype=np.double,order='F')
            double[:,:] uflxc_lw_out = np.zeros((ncol,nz + n_ext +1),dtype=np.double,order='F')
            double[:,:] dflxc_lw_out = np.zeros((ncol,nz + n_ext +1),dtype=np.double,order='F')
            double[:,:] hrc_lw_out = np.zeros((ncol,nz + n_ext),dtype=np.double,order='F')
            double[:,:] duflx_dt_out = np.zeros((ncol,nz + n_ext +1),dtype=np.double,order='F')
            double[:,:] duflxc_dt_out = np.zeros((ncol,nz + n_ext +1),dtype=np.double,order='F')
            double[:,:] uflx_sw_out = np.zeros((ncol,nz + n_ext +1),dtype=np.double,order='F')
            double[:,:] dflx_sw_out = np.zeros((ncol,nz + n_ext +1),dtype=np.double,order='F')
            double[:,:] hr_sw_out = np.zeros((ncol,nz + n_ext),dtype=np.double,order='F')
            double[:,:] uflxc_sw_out = np.zeros((ncol,nz + n_ext +1),dtype=np.double,order='F')
            double[:,:] dflxc_sw_out = np.zeros((ncol,nz + n_ext +1),dtype=np.double,order='F')
            double[:,:] hrc_sw_out = np.zeros((ncol,nz + n_ext),dtype=np.double,order='F')
            
        
        # toc=timeit.default_timer()
        # print 'Timer 4: ', toc-tic
        # tic=toc
        cdef double rv_to_reff = np.exp(np.log(1.2)**2.0)*10.0*1000.0
        
        with nogil:
            for i in prange(nxl-2*gw,schedule='static'):
            # for i in xrange(nxl-2*gw):
                for j in xrange(nyl-2*gw):
                    ishift = i + gw
                    jshift = j + gw
                    ijshift = i*(nyl-2*gw)+j
                    
                    tsfc_in[ijshift] = surface_temperature
                    coszen_in[ijshift] = coszen
                    asdif_in[ijshift]  = adif
                    aldif_in[ijshift]  = adif
                    asdir_in[ijshift]  = adir
                    aldir_in[ijshift]  = adir
                    for k in xrange(16):
                        emis_in[ijshift, k] = 0.95
                    
                    for k in xrange(nz + n_ext):
                        kshift = k + gw
                        
                        play_in[ijshift, k] = pressure_ext[k]/100.0
                        plev_in[ijshift, k] = p0i_ext[k]/100.0
                        tlay_in[ijshift, k] = temperature_ext[i, j, k]
                        # tlev_in: a separate loop
                        h2ovmr_in[ijshift, k] = Rv/Ra*y_vapor_ext[i,j,k] * h2o_factor # Irrational Test Parameter...
                        o3vmr_in [ijshift, k] = o3vmr[k]
                        co2vmr_in[ijshift, k] = co2vmr[k]
                        ch4vmr_in[ijshift, k] = ch4vmr[k]
                        n2ovmr_in[ijshift, k] = n2ovmr[k]
                        o2vmr_in [ijshift, k] = o2vmr[k]
                        o3vmr_in [ijshift, k] = o3vmr[k]
                        cfc11vmr_in[ijshift, k] = cfc11vmr[k]
                        cfc12vmr_in[ijshift, k] = cfc12vmr[k]
                        cfc22vmr_in[ijshift, k] = cfc22vmr[k]
                        ccl4vmr_in[ijshift, k] = ccl4vmr[k]
                        cldfr_in[ijshift, k] = cloud_fraction_ext[i,j,k]
                        cliqwp_in[ijshift, k] = y_cond_ext[i,j,k]*1.0e3*(p0i_ext[k] - p0i_ext[k+1])/g
                        if uniform_reliq:
                            reliq_in[ijshift, k] = 14.0*cloud_fraction_ext[i,j,k]
                        else:
                            reliq_in[ijshift, k] = ((3.0*pressure_ext[k]/Ra/temperature_ext[i,j,k]*y_cond_ext[i,j,k]/
                                                    fmax(cloud_fraction_ext[i,j,k],1.0e-6))/(4.0*pi*1.0e3*100.0))**(1.0/3.0)
                            reliq_in[ijshift, k] = fmin(fmax(reliq_in[ijshift, k]*rv_to_reff, 2.5), 60.0)
                            
                    tlev_in[ijshift, 0] = surface_temperature
                    for k in xrange(1,nz+n_ext):
                        tlev_in[ijshift, k] = 0.5*(tlay_in[ijshift,k-1]+tlay_in[ijshift,k])
                    tlev_in[ijshift, nz+n_ext] = 2.0*tlay_in[ijshift,nz+n_ext-1] - tlev_in[ijshift,nz+n_ext-1]
        
        # toc=timeit.default_timer()
        # print 'Timer 5: ', toc-tic
        # tic=toc
        
        # BEGIN: Debug Outputs 
        cdef double [:] play_in_ext_   = np.zeros(nz+ n_ext,dtype=np.double,order='F')
        cdef double [:] plev_in_ext_   = np.zeros(nz+ n_ext +1,dtype=np.double,order='F')
        cdef double [:] tlay_in_ext_   = np.zeros(nz+ n_ext,dtype=np.double,order='F')
        cdef double [:] tlev_in_ext_   = np.zeros(nz+ n_ext +1,dtype=np.double,order='F')
        cdef double [:] h2ovmr_in_ext_ = np.zeros(nz+ n_ext,dtype=np.double,order='F')
        cdef double tsfc_in_ext_ = 0.0
        cdef double [:] cldfr_in_ext_  = np.zeros(nz+ n_ext,dtype=np.double,order='F')
        cdef double [:] cliqwp_in_ext_ = np.zeros(nz+ n_ext,dtype=np.double,order='F')
        cdef double [:] reliq_in_ext_  = np.zeros(nz+ n_ext,dtype=np.double,order='F')
        
        with nogil:
            for i in xrange(ncol): 
            # ncol = (nxl-2*gw)*(nyl-2*gw)
                for k in xrange(nz+n_ext+1):
                    plev_in_ext_[k] += plev_in[i,k] # Debug Outputs
                    tlev_in_ext_[k] += tlev_in[i,k] # Debug Outputs
                    
                for k in xrange(nz+n_ext):
                    play_in_ext_[k] += play_in[i,k] # Debug Outputs
                    tlay_in_ext_[k] += tlay_in[i,k] # Debug Outputs
                    h2ovmr_in_ext_[k] += h2ovmr_in[i,k] # Debug Outputs
                    cldfr_in_ext_[k]  += cldfr_in[i,k] # Debug Outputs
                    cliqwp_in_ext_[k] += cliqwp_in[i,k] # Debug Outputs
                    reliq_in_ext_[k]  += reliq_in[i,k] # Debug Outputs
                    
                tsfc_in_ext_ += tsfc_in[i]
        
        cdef double [:] play_in_ext   = np.zeros(nz+ n_ext,dtype=np.double,order='F')
        cdef double [:] plev_in_ext   = np.zeros(nz+ n_ext +1,dtype=np.double,order='F')
        cdef double [:] tlay_in_ext   = np.zeros(nz+ n_ext,dtype=np.double,order='F')
        cdef double [:] tlev_in_ext   = np.zeros(nz+ n_ext +1,dtype=np.double,order='F')
        cdef double [:] h2ovmr_in_ext = np.zeros(nz+ n_ext,dtype=np.double,order='F')
        cdef double tsfc_in_ext = 0.0
        cdef double [:] cldfr_in_ext  = np.zeros(nz+ n_ext,dtype=np.double,order='F')
        cdef double [:] cliqwp_in_ext = np.zeros(nz+ n_ext,dtype=np.double,order='F')
        cdef double [:] reliq_in_ext  = np.zeros(nz+ n_ext,dtype=np.double,order='F')
        
        comm.cart_comm.Allreduce(play_in_ext_,   play_in_ext,    op = MPI.SUM)
        comm.cart_comm.Allreduce(plev_in_ext_,   plev_in_ext,    op = MPI.SUM)
        comm.cart_comm.Allreduce(tlay_in_ext_,   tlay_in_ext,    op = MPI.SUM)
        comm.cart_comm.Allreduce(tlev_in_ext_,   tlev_in_ext,    op = MPI.SUM)
        comm.cart_comm.Allreduce(h2ovmr_in_ext_, h2ovmr_in_ext,  op = MPI.SUM)
        comm.cart_comm.Allreduce(cldfr_in_ext_,  cldfr_in_ext,   op = MPI.SUM)
        comm.cart_comm.Allreduce(cliqwp_in_ext_, cliqwp_in_ext,  op = MPI.SUM)
        comm.cart_comm.Allreduce(reliq_in_ext_,  reliq_in_ext,   op = MPI.SUM)
        #comm.cart_comm.Allreduce(tsfc_in_ext_,   tsfc_in_ext,    op = MPI.SUM)    
            
        with nogil:
            for k in xrange(nz+n_ext+1):
            # nxny = (nx * ny)
                plev_in_ext[k] = plev_in_ext[k]/nxny # Debug Outputs
                tlev_in_ext[k] = tlev_in_ext[k]/nxny # Debug Outputs
            for k in xrange(nz+n_ext):
                play_in_ext[k]  = play_in_ext[k]/nxny    # Debug Outputs
                tlay_in_ext[k]  = tlay_in_ext[k]/nxny    # Debug Outputs
                h2ovmr_in_ext[k] = h2ovmr_in_ext[k]/nxny # Debug Outputs
                cldfr_in_ext[k]  = cldfr_in_ext[k]/nxny  # Debug Outputs
                cliqwp_in_ext[k] = cliqwp_in_ext[k]/nxny # Debug Outputs
                reliq_in_ext[k]  = reliq_in_ext[k]/nxny # Debug Outputs
            tsfc_in_ext = surface_temperature #tsfc_in_ext/nxny
            
        if self.debug_mode:
            self.play_in_ext  = np.copy(play_in_ext)
            self.plev_in_ext  = np.copy(plev_in_ext)      
            self.tlay_in_ext  = np.copy(tlay_in_ext)     
            self.tlev_in_ext  = np.copy(tlev_in_ext)     
            self.h2ovmr_in_ext  = np.copy(h2ovmr_in_ext)
            self.cldfr_in_ext   = np.copy(cldfr_in_ext)
            self.cliqwp_in_ext  = np.copy(cliqwp_in_ext)
            self.reliq_in_ext   = np.copy(reliq_in_ext)
            self.tsfc_in_ext   = np.copy(tsfc_in_ext)
            
            if comm.rank == 0:
                tic=timeit.default_timer()
        # toc=timeit.default_timer()
        # print 'Timer 6: ', toc-tic
        # tic=toc
            
        # np.set_printoptions(threshold=np.nan)
        # if comm.rank == 0:
            # print 'p0i_ext',       self.p0i_ext[:]
            # print 'pressure_ext',  self.pressure_ext[:]
            # print 'play_in_ext',   self.play_in_ext[:]
            # print 'plev_in_ext',   self.plev_in_ext[:]
            # print 'tlay_in_ext',   self.tlay_in_ext[:]
            # print 'tlev_in_ext',   self.tlev_in_ext[:]
            # print 'h2ovmr_in_ext', self.h2ovmr_in_ext[:]
            # print 'cldfr_in_ext',  self.cldfr_in_ext[:]
            # print 'cliqwp_in_ext', self.cliqwp_in_ext[:]
            # print 'tsfc_in_ext',   self.tsfc_in_ext
        # END: Debug Outputs 
        
        c_rrtmg_lw (
             &ncol    ,&nlay    ,&icld    ,&idrv    , 
             &play_in[0,0]    ,&plev_in[0,0]    ,&tlay_in[0,0]    ,&tlev_in[0,0]    ,&tsfc_in[0]    , 
             &h2ovmr_in[0,0]  ,&o3vmr_in[0,0]   ,&co2vmr_in[0,0]  ,&ch4vmr_in[0,0]  ,&n2ovmr_in[0,0]  ,&o2vmr_in[0,0], 
             &cfc11vmr_in[0,0],&cfc12vmr_in[0,0],&cfc22vmr_in[0,0],&ccl4vmr_in[0,0] ,&emis_in[0,0]    , 
             &inflglw ,&iceflglw,&liqflglw,&cldfr_in[0,0]   , 
             &taucld_lw_in[0,0,0]  ,&cicewp_in[0,0]  ,&cliqwp_in[0,0]  ,&reice_in[0,0]   ,&reliq_in[0,0]   , 
             &tauaer_lw_in[0,0,0]  , 
             &uflx_lw_out[0,0]    ,&dflx_lw_out[0,0]    ,&hr_lw_out[0,0]      ,&uflxc_lw_out[0,0]   ,&dflxc_lw_out[0,0],  &hrc_lw_out[0,0], 
             &duflx_dt_out[0,0],&duflxc_dt_out[0,0] )        

        
        # toc=timeit.default_timer()
        # print 'Timer 7 (RRTM LW): ', toc-tic
        # tic=toc
        
        c_rrtmg_sw (
             &ncol    ,&nlay    ,&icld    ,&iaer    , 
             &play_in[0,0]    ,&plev_in[0,0]    ,&tlay_in[0,0]    ,&tlev_in[0,0]    ,&tsfc_in[0]    , 
             &h2ovmr_in[0,0]  ,&o3vmr_in[0,0]   ,&co2vmr_in[0,0]  ,&ch4vmr_in[0,0]  ,&n2ovmr_in[0,0]  ,&o2vmr_in[0,0], 
             &asdir_in[0]   ,&asdif_in[0]   ,&aldir_in[0]   ,&aldif_in[0]   , 
             &coszen_in[0]  ,&adjes_in   ,&dyofyr_in  ,&scon_in    , 
             &inflgsw ,&iceflgsw,&liqflgsw,&cldfr_in[0,0]   , 
             &taucld_sw_in[0,0,0]  ,&ssacld_sw_in[0,0,0]  ,&asmcld_sw_in[0,0,0]  ,&fsfcld_sw_in[0,0,0]  , 
             &cicewp_in[0,0]  ,&cliqwp_in[0,0]  ,&reice_in[0,0]   ,&reliq_in[0,0]   , 
             &tauaer_sw_in[0,0,0]  ,&ssaaer_sw_in[0,0,0]  ,&asmaer_sw_in[0,0,0]  ,&ecaer_sw_in[0,0,0]   , 
             &uflx_sw_out[0,0]    ,&dflx_sw_out[0,0]    ,&hr_sw_out[0,0]      ,&uflxc_sw_out[0,0]   ,&dflxc_sw_out[0,0],  &hrc_sw_out[0,0])
             
        
        # toc=timeit.default_timer()
        # print 'Timer 8 (RRTM SW): ', toc-tic
        # tic=toc
        
        cdef double [:,:] uflx_lw  = np.asarray(uflx_lw_out, dtype=np.double,order='F')
        cdef double [:,:] dflx_lw  = np.asarray(dflx_lw_out, dtype=np.double,order='F')
        cdef double [:,:] hr_lw    = np.asarray(hr_lw_out,   dtype=np.double,order='F')
        cdef double [:,:] uflxc_lw = np.asarray(uflxc_lw_out,dtype=np.double,order='F')
        cdef double [:,:] dflxc_lw = np.asarray(dflxc_lw_out,dtype=np.double,order='F')
        cdef double [:,:] hrc_lw   = np.asarray(hrc_lw_out,  dtype=np.double,order='F')
        cdef double [:,:] uflx_sw  = np.asarray(uflx_sw_out, dtype=np.double,order='F') 
        cdef double [:,:] dflx_sw  = np.asarray(dflx_sw_out, dtype=np.double,order='F') 
        cdef double [:,:] hr_sw    = np.asarray(hr_sw_out,   dtype=np.double,order='F') 
        cdef double [:,:] uflxc_sw = np.asarray(uflxc_sw_out,dtype=np.double,order='F') 
        cdef double [:,:] dflxc_sw = np.asarray(dflxc_sw_out,dtype=np.double,order='F') 
        cdef double [:,:] hrc_sw   = np.asarray(hrc_sw_out,  dtype=np.double,order='F')
        
        # Averaging using MPI ...
        # Flux output dim: (ncol,nz + n_ext +1); heat-rate output dim: (ncol,nz + n_ext)
        
        cdef double [:] lw_up_flux_ext_    = np.zeros(nz+ n_ext+1,dtype=np.double,order='F')
        cdef double [:] lw_down_flux_ext_  = np.zeros(nz+ n_ext+1,dtype=np.double,order='F')
        cdef double [:] sw_up_flux_ext_    = np.zeros(nz+ n_ext+1,dtype=np.double,order='F')
        cdef double [:] sw_down_flux_ext_  = np.zeros(nz+ n_ext+1,dtype=np.double,order='F')
        cdef double [:] lw_up_fluxc_ext_   = np.zeros(nz+ n_ext+1,dtype=np.double,order='F')
        cdef double [:] lw_down_fluxc_ext_ = np.zeros(nz+ n_ext+1,dtype=np.double,order='F')
        cdef double [:] sw_up_fluxc_ext_   = np.zeros(nz+ n_ext+1,dtype=np.double,order='F')
        cdef double [:] sw_down_fluxc_ext_ = np.zeros(nz+ n_ext+1,dtype=np.double,order='F')
        cdef double [:] lw_tt_ext_    = np.zeros(nz+ n_ext,dtype=np.double,order='F')
        cdef double [:] sw_tt_ext_    = np.zeros(nz+ n_ext,dtype=np.double,order='F')
        cdef double [:] total_tt_ext_ = np.zeros(nz+ n_ext,dtype=np.double,order='F')
        cdef double [:] lw_ttc_ext_    = np.zeros(nz+ n_ext,dtype=np.double,order='F')
        cdef double [:] sw_ttc_ext_    = np.zeros(nz+ n_ext,dtype=np.double,order='F')
        cdef double [:] total_ttc_ext_ = np.zeros(nz+ n_ext,dtype=np.double,order='F')
        
        
        with nogil:
            for i in xrange(ncol): 
            # ncol = (nxl-2*gw)*(nyl-2*gw)
                for k in xrange(nz+n_ext+1):
                    lw_up_flux_ext_[k]    +=  uflx_lw[i,k]
                    lw_down_flux_ext_[k]  +=  dflx_lw[i,k]
                    sw_up_flux_ext_[k]    +=  uflx_sw[i,k]
                    sw_down_flux_ext_[k]  +=  dflx_sw[i,k]
                    lw_up_fluxc_ext_[k]   +=  uflxc_lw[i,k]
                    lw_down_fluxc_ext_[k] +=  dflxc_lw[i,k]
                    sw_up_fluxc_ext_[k]   +=  uflxc_sw[i,k]
                    sw_down_fluxc_ext_[k] +=  dflxc_sw[i,k]
                    
                for k in xrange(nz+n_ext):
                    lw_tt_ext_[k]  +=  hr_lw[i,k]
                    sw_tt_ext_[k]  +=  hr_sw[i,k]
                    lw_ttc_ext_[k] +=  hrc_lw[i,k]
                    sw_ttc_ext_[k] +=  hrc_sw[i,k]
                    total_tt_ext_[k]  += (hr_lw[i,k]+hr_sw[i,k])
                    total_ttc_ext_[k] += (hrc_lw[i,k]+hrc_sw[i,k])
                
        
        cdef double [:] lw_up_flux_ext    = np.zeros(nz+ n_ext+1,dtype=np.double,order='F')
        cdef double [:] lw_down_flux_ext  = np.zeros(nz+ n_ext+1,dtype=np.double,order='F')
        cdef double [:] sw_up_flux_ext    = np.zeros(nz+ n_ext+1,dtype=np.double,order='F')
        cdef double [:] sw_down_flux_ext  = np.zeros(nz+ n_ext+1,dtype=np.double,order='F')
        cdef double [:] lw_up_fluxc_ext   = np.zeros(nz+ n_ext+1,dtype=np.double,order='F')
        cdef double [:] lw_down_fluxc_ext = np.zeros(nz+ n_ext+1,dtype=np.double,order='F')
        cdef double [:] sw_up_fluxc_ext   = np.zeros(nz+ n_ext+1,dtype=np.double,order='F')
        cdef double [:] sw_down_fluxc_ext = np.zeros(nz+ n_ext+1,dtype=np.double,order='F')
        cdef double [:] lw_tt_ext    = np.zeros(nz+ n_ext,dtype=np.double,order='F')
        cdef double [:] sw_tt_ext    = np.zeros(nz+ n_ext,dtype=np.double,order='F')
        cdef double [:] total_tt_ext = np.zeros(nz+ n_ext,dtype=np.double,order='F')
        cdef double [:] lw_ttc_ext    = np.zeros(nz+ n_ext,dtype=np.double,order='F')
        cdef double [:] sw_ttc_ext    = np.zeros(nz+ n_ext,dtype=np.double,order='F')
        cdef double [:] total_ttc_ext = np.zeros(nz+ n_ext,dtype=np.double,order='F')
        
        comm.cart_comm.Allreduce(lw_up_flux_ext_,   lw_up_flux_ext,  op = MPI.SUM)
        comm.cart_comm.Allreduce(lw_down_flux_ext_, lw_down_flux_ext,op = MPI.SUM)
        comm.cart_comm.Allreduce(sw_up_flux_ext_,   sw_up_flux_ext,  op = MPI.SUM)
        comm.cart_comm.Allreduce(sw_down_flux_ext_, sw_down_flux_ext,op = MPI.SUM)
        comm.cart_comm.Allreduce(lw_up_fluxc_ext_,   lw_up_fluxc_ext,  op = MPI.SUM)
        comm.cart_comm.Allreduce(lw_down_fluxc_ext_, lw_down_fluxc_ext,op = MPI.SUM)
        comm.cart_comm.Allreduce(sw_up_fluxc_ext_,   sw_up_fluxc_ext,  op = MPI.SUM)
        comm.cart_comm.Allreduce(sw_down_fluxc_ext_, sw_down_fluxc_ext,op = MPI.SUM)
        comm.cart_comm.Allreduce(lw_tt_ext_,    lw_tt_ext,    op = MPI.SUM)
        comm.cart_comm.Allreduce(sw_tt_ext_,    sw_tt_ext,    op = MPI.SUM)
        comm.cart_comm.Allreduce(total_tt_ext_, total_tt_ext, op = MPI.SUM)
        comm.cart_comm.Allreduce(lw_ttc_ext_,    lw_ttc_ext,    op = MPI.SUM)
        comm.cart_comm.Allreduce(sw_ttc_ext_,    sw_ttc_ext,    op = MPI.SUM)
        comm.cart_comm.Allreduce(total_ttc_ext_, total_ttc_ext, op = MPI.SUM)
        
        with nogil:
            for k in xrange(nz+n_ext+1):
            # nxny = (nx * ny)
                lw_up_flux_ext[k]    = lw_up_flux_ext[k]   /nxny
                lw_down_flux_ext[k]  = lw_down_flux_ext[k] /nxny
                sw_up_flux_ext[k]    = sw_up_flux_ext[k]   /nxny
                sw_down_flux_ext[k]  = sw_down_flux_ext[k] /nxny
                lw_up_fluxc_ext[k]   = lw_up_fluxc_ext[k]  /nxny
                lw_down_fluxc_ext[k] = lw_down_fluxc_ext[k]/nxny
                sw_up_fluxc_ext[k]   = sw_up_fluxc_ext[k]  /nxny
                sw_down_fluxc_ext[k] = sw_down_fluxc_ext[k]/nxny
            for k in xrange(nz+n_ext):
                # Convert T-tendencies from K/day to K/s
                lw_tt_ext[k]     = lw_tt_ext[k]    /nxny /86400.0
                sw_tt_ext[k]     = sw_tt_ext[k]    /nxny /86400.0
                total_tt_ext[k]  = total_tt_ext[k] /nxny /86400.0
                lw_ttc_ext[k]    = lw_ttc_ext[k]   /nxny /86400.0
                sw_ttc_ext[k]    = sw_ttc_ext[k]   /nxny /86400.0
                total_ttc_ext[k] = total_ttc_ext[k]/nxny /86400.0
   
        
        self.lw_flux_down = np.copy(lw_down_flux_ext[:nz])
        self.lw_flux_up   = np.copy(lw_up_flux_ext[:nz])
        self.sw_flux_down = np.copy(sw_down_flux_ext[:nz])
        self.sw_flux_up   = np.copy(sw_up_flux_ext[:nz])

        self.lw_flux_down_clear = np.copy(lw_down_fluxc_ext[:nz])
        self.lw_flux_up_clear   = np.copy(lw_up_fluxc_ext[:nz])
        self.sw_flux_down_clear = np.copy(sw_down_fluxc_ext[:nz])
        self.sw_flux_up_clear   = np.copy(sw_up_fluxc_ext[:nz])
        self.dTdt_clear = np.copy(total_ttc_ext[:nz])
        
        cdef double [:,:] osr  = self.osr
        cdef double [:,:] olr  = self.olr
        with nogil:
            for i in prange(nxl-2*gw,schedule='static'):
                for j in xrange(nyl-2*gw):
                    ijshift = i*(nyl-2*gw)+j
                    osr[i+gw,j+gw] = uflx_sw[ijshift,nz]
                    olr[i+gw,j+gw] = uflx_lw[ijshift,nz]
                
        self.total_tt_ext = np.copy(total_tt_ext)
        
        # Added for 3D output
        with nogil:
            for i in prange(nxl-2*gw,schedule='static'):
            # for i in xrange(nxl-2*gw):
                for j in xrange(nyl-2*gw):
                    ijshift = i*(nyl-2*gw)+j
                    for k in xrange(nz + n_ext):
                        total_tt_ext_3d[i, j, k] = (hr_lw[ijshift,k]+hr_sw[ijshift,k])/86400.0
                        
        self.total_tt_ext_3d = np.copy(total_tt_ext_3d)
        
                    
                    
        # DEBUG OUTPUTS 
        if self.debug_mode:
            self.lw_up_flux_ext   = np.copy(lw_up_flux_ext)
            self.lw_down_flux_ext = np.copy(lw_down_flux_ext)
            self.sw_up_flux_ext   = np.copy(sw_up_flux_ext)
            self.sw_down_flux_ext = np.copy(sw_down_flux_ext)
            self.lw_up_fluxc_ext   = np.copy(lw_up_fluxc_ext)
            self.lw_down_fluxc_ext = np.copy(lw_down_fluxc_ext)
            self.sw_up_fluxc_ext   = np.copy(sw_up_fluxc_ext)
            self.sw_down_fluxc_ext = np.copy(sw_down_fluxc_ext)
            self.lw_tt_ext     = np.copy(lw_tt_ext)
            self.sw_tt_ext     = np.copy(sw_tt_ext)
            self.lw_ttc_ext    = np.copy(lw_ttc_ext)
            self.sw_ttc_ext    = np.copy(sw_ttc_ext)
            self.total_ttc_ext = np.copy(total_ttc_ext)
            
        # if comm.rank == 0:
            # print 'lw_up_flux_ext',     self.lw_up_flux_ext[:]
            # print 'lw_down_flux_ext',   self.lw_down_flux_ext[:]
            # print 'sw_up_flux_ext',     self.sw_up_flux_ext[:]
            # print 'sw_down_flux_ext',   self.sw_down_flux_ext[:]
            # print 'lw_up_fluxc_ext',    self.lw_up_fluxc_ext[:]
            # print 'lw_down_fluxc_ext',  self.lw_down_fluxc_ext[:]
            # print 'sw_up_fluxc_ext',    self.sw_up_fluxc_ext[:]
            # print 'sw_down_fluxc_ext',  self.sw_down_fluxc_ext[:]
            # print 'lw_tt_ext',          self.lw_tt_ext[:]
            # print 'sw_tt_ext',          self.sw_tt_ext[:]
            #print 'total_tt_ext',       self.total_tt_ext[:]
            # print 'lw_ttc_ext',         self.lw_ttc_ext[:]
            # print 'sw_ttc_ext',         self.sw_ttc_ext[:]
            # print 'total_ttc_ext',      self.total_ttc_ext[:]
            
            if comm.rank == 0:
                toc=timeit.default_timer()
                print 'RRTM Computing Time: ', toc-tic
        # tic=toc
        
        
        return

    def init_output(self,grid,io):
        io.init_profile('lw_flux_down','wm^-2',grid.nz)
        io.init_profile('lw_flux_up','wm^-2',grid.nz)
        io.init_profile('sw_flux_down','wm^-2',grid.nz)
        io.init_profile('sw_flux_up','wm^-2',grid.nz)
        io.init_profile('dTdt_rad','Ks^-1',grid.nz)
        io.init_profile('dsdt_rad','J(kgK)^-1s^-1',grid.nz)
        # New outputs:
        io.init_profile('lw_flux_down_clear','wm^-2',grid.nz)
        io.init_profile('lw_flux_up_clear','wm^-2',grid.nz)
        io.init_profile('sw_flux_down_clear','wm^-2',grid.nz)
        io.init_profile('sw_flux_up_clear','wm^-2',grid.nz)
        io.init_profile('dTdt_rad_clear','Ks^-1',grid.nz)
        
        if io.output_2d_field:
            io.init_topview('osr', 'wm^-2', grid.nx, grid.ny)
            io.init_topview('olr', 'wm^-2', grid.nx, grid.ny)
            
            
        # Other radiation parameters: output at initial step only
        io.init_timeseries('dyofyr','-')
        io.init_timeseries('adjes','-')
        io.init_timeseries('scon','-')
        io.init_timeseries('coszen','-')
        io.init_timeseries('adif','-')
        io.init_timeseries('adir','-')  
           
        try:
            io.output_timeseries(self.dyofyr,'dyofyr')
        except:
            print('dyofyr error') 
        try:
            io.output_timeseries(self.adjes,'adjes')
        except:
            print('adjes error')  
        try:
            io.output_timeseries(self.scon,'scon')
        except:
            print('scon error')  
        try:
            io.output_timeseries(self.coszen,'coszen')
        except:
            print('coszen error')  
        try:
            io.output_timeseries(self.adif,'adif')
        except:
            print('adif error')  
        try:
            io.output_timeseries(self.adir,'adir')
        except:
            print('adir error')  
                          
        # DEBUG OUTPUTS INITIALIZATION
        # print (grid.nz+self.n_ext)
        if self.debug_mode:
            io.init_profile('p0i_ext','Pa',grid.nz+self.n_ext+1)
            io.init_profile('pressure_ext','Pa',grid.nz+self.n_ext)
            io.init_profile('play_in_ext','hPa',grid.nz+self.n_ext)
            io.init_profile('plev_in_ext','hPa',grid.nz+self.n_ext+1)
            io.init_profile('tlay_in_ext','K',grid.nz+self.n_ext)
            io.init_profile('tlev_in_ext','K',grid.nz+self.n_ext+1)
            io.init_profile('h2ovmr_in_ext','1',grid.nz+self.n_ext)
            io.init_profile('cldfr_in_ext', '1',grid.nz+self.n_ext)
            io.init_profile('cliqwp_in_ext','gm^-2',grid.nz+self.n_ext)
            io.init_profile('reliq_in_ext','micron',grid.nz+self.n_ext)
            io.init_profile('lw_up_flux_ext',  'wm^-2',grid.nz+self.n_ext+1)
            io.init_profile('lw_down_flux_ext','wm^-2',grid.nz+self.n_ext+1)
            io.init_profile('sw_up_flux_ext',  'wm^-2',grid.nz+self.n_ext+1)
            io.init_profile('sw_down_flux_ext','wm^-2',grid.nz+self.n_ext+1)
            io.init_profile('lw_up_fluxc_ext',  'wm^-2',grid.nz+self.n_ext+1)
            io.init_profile('lw_down_fluxc_ext','wm^-2',grid.nz+self.n_ext+1)
            io.init_profile('sw_up_fluxc_ext',  'wm^-2',grid.nz+self.n_ext+1)
            io.init_profile('sw_down_fluxc_ext','wm^-2',grid.nz+self.n_ext+1)
            io.init_profile('lw_tt_ext',   'K/s',grid.nz+self.n_ext)
            io.init_profile('sw_tt_ext',   'K/s',grid.nz+self.n_ext)
            io.init_profile('total_tt_ext','K/s',grid.nz+self.n_ext)
            io.init_profile('lw_ttc_ext',   'K/s',grid.nz+self.n_ext)
            io.init_profile('sw_ttc_ext',   'K/s',grid.nz+self.n_ext)
            io.init_profile('total_ttc_ext','K/s',grid.nz+self.n_ext)
            io.init_timeseries('tsfc_in_ext','K')
            # Trace Gases: output at initial step only
            io.init_profile('o3vmr','-',grid.nz+self.n_ext)
            io.init_profile('co2vmr','-',grid.nz+self.n_ext)
            io.init_profile('ch4vmr','-',grid.nz+self.n_ext)
            io.init_profile('n2ovmr','-',grid.nz+self.n_ext)
            io.init_profile('o2vmr','-',grid.nz+self.n_ext)
            io.init_profile('cfc11vmr','-',grid.nz+self.n_ext)
            io.init_profile('cfc12vmr','-',grid.nz+self.n_ext)
            io.init_profile('cfc22vmr','-',grid.nz+self.n_ext)
            io.init_profile('ccl4vmr','-',grid.nz+self.n_ext)

            
            #Output extended T/p profiles
            try:
                io.output_profile(self.pressure_ext,'pressure_ext')
            except:
                print('pressure_ext error')    
            try:
                io.output_profile(self.p0i_ext,'p0i_ext')
            except:
                print('p0i_ext error')    
            try:
                io.output_profile(self.o3vmr,'o3vmr')
            except:
                print('o3vmr error')     
            try:
                io.output_profile(self.co2vmr,'co2vmr')
            except:
                print('co2vmr error')     
            try:
                io.output_profile(self.ch4vmr,'ch4vmr')
            except:
                print('ch4vmr error')     
            try:
                io.output_profile(self.n2ovmr,'n2ovmr')
            except:
                print('n2ovmr error')     
            try:
                io.output_profile(self.o2vmr,'o2vmr')
            except:
                print('o2vmr error')     
            try:
                io.output_profile(self.cfc11vmr,'cfc11vmr')
            except:
                print('cfc11vmr error') 
            try:
                io.output_profile(self.cfc12vmr,'cfc12vmr')
            except:
                print('cfc12vmr error') 
            try:
                io.output_profile(self.cfc22vmr,'cfc22vmr')
            except:
                print('cfc22vmr error') 
            try:
                io.output_profile(self.ccl4vmr,'ccl4vmr')
            except:
                print('ccl4vmr error') 

                                
        return

    def output(self,grid,io,comm):
    
        #Output lw fluxes
        try:
            io.output_profile(self.lw_flux_down,'lw_flux_down')
        except:
            if comm.rank == 0:
                print('lw_flux_down error')
        try:
            io.output_profile(self.lw_flux_up,'lw_flux_up')
        except:
            if comm.rank == 0:
                print('lw_flux_up error')
            
        #Output sw fluxes
        try:
            io.output_profile(self.sw_flux_down,'sw_flux_down')
        except:
            if comm.rank == 0:
                print('sw_flux_down error')
        try:
            io.output_profile(self.sw_flux_up,'sw_flux_up')
        except:
            if comm.rank == 0:
                print('sw_flux_up error')
            
            
        #Output T and s tendencies
        try:
            _mean, _mean_squared, _variance, _mean_cubed, _skew  =   ProfileStatistics.resolved_statistics(self.dTdt_rad,grid,comm)
            io.output_profile(_mean, 'dTdt_rad')
        except:
            if comm.rank == 0:
                print('dTdt_rad error')
            
        try:
            _mean, _mean_squared, _variance, _mean_cubed, _skew  =   ProfileStatistics.resolved_statistics(self.dsdt_rad,grid,comm)
            io.output_profile(_mean, 'dsdt_rad')
        except:
            if comm.rank == 0:
                print('dsdt_rad error')
            
        #Clear Outputs
        #Output lw fluxes
        try:
            io.output_profile(self.lw_flux_down_clear,'lw_flux_down_clear')
        except:
            if comm.rank == 0:
                print('lw_flux_down_clear error')
        try:
            io.output_profile(self.lw_flux_up_clear,'lw_flux_up_clear')
        except:
            if comm.rank == 0:
                print('lw_flux_up_clear error')
            
        #Output sw fluxes
        try:
            io.output_profile(self.sw_flux_down_clear,'sw_flux_down_clear')
        except:
            if comm.rank == 0:
                print('sw_flux_down_clear error')
        try:
            io.output_profile(self.sw_flux_up_clear,'sw_flux_up_clear')
        except:
            if comm.rank == 0:
                print('sw_flux_up_clear error')
        try:
            io.output_profile(self.dTdt_clear,'dTdt_rad_clear')
        except:
            if comm.rank == 0:
                print('dTdt_rad_clear error')
        
        if io.output_2d_field:
            reduced_data = reduce_to_root_2d(self.osr,grid,comm)
            io.output_topview(reduced_data,'osr')
            reduced_data = reduce_to_root_2d(self.olr,grid,comm)
            io.output_topview(reduced_data,'olr')
        	# Output into field file
        	
                            
        # DEBUG OUTPUTS BELOW
        if self.debug_mode:
   
            try:
                io.output_profile(self.play_in_ext,'play_in_ext')
            except:
                if comm.rank == 0:
                    print('play_in_ext error')    
            try:
                io.output_profile(self.plev_in_ext,'plev_in_ext')
            except:
                if comm.rank == 0:
                    print('plev_in_ext error')    
            try:
                io.output_profile(self.tlay_in_ext,'tlay_in_ext')
            except:
                if comm.rank == 0:
                    print('tlay_in_ext error')    
            try:
                io.output_profile(self.tlev_in_ext,'tlev_in_ext')
            except:
                if comm.rank == 0:
                    print('tlev_in_ext error')    
            try:
                io.output_profile(self.h2ovmr_in_ext,'h2ovmr_in_ext')
            except:
                if comm.rank == 0:
                    print('h2ovmr_in_ext error')    
            try:
                io.output_profile(self.cldfr_in_ext,'cldfr_in_ext')
            except:
                if comm.rank == 0:
                    print('cldfr_in_ext error')    
            try:
                io.output_profile(self.cliqwp_in_ext,'cliqwp_in_ext')
            except:
                if comm.rank == 0:
                    print('cliqwp_in_ext error')    
            try:
                io.output_profile(self.reliq_in_ext,'reliq_in_ext')
            except:
                if comm.rank == 0:
                    print('reliq_in_ext error') 
            try:
                io.output_profile(self.lw_up_flux_ext,'lw_up_flux_ext')
            except:
                if comm.rank == 0:
                    print('lw_up_flux_ext error')    
            try:
                io.output_profile(self.lw_down_flux_ext,'lw_down_flux_ext')
            except:
                if comm.rank == 0:
                    print('lw_down_flux_ext error')    
            try:
                io.output_profile(self.sw_up_flux_ext,'sw_up_flux_ext')
            except:
                if comm.rank == 0:
                    print('sw_up_flux_ext error')    
            try:
                io.output_profile(self.sw_down_flux_ext,'sw_down_flux_ext')
            except:
                if comm.rank == 0:
                    print('sw_down_flux_ext error')    
            try:
                io.output_profile(self.lw_up_fluxc_ext,'lw_up_fluxc_ext')
            except:
                if comm.rank == 0:
                    print('lw_up_fluxc_ext error')    
            try:
                io.output_profile(self.lw_down_fluxc_ext,'lw_down_fluxc_ext')
            except:
                if comm.rank == 0:
                    print('lw_down_fluxc_ext error')    
            try:
                io.output_profile(self.sw_up_fluxc_ext,'sw_up_fluxc_ext')
            except:
                if comm.rank == 0:
                    print('sw_up_fluxc_ext error')    
            try:
                io.output_profile(self.sw_down_fluxc_ext,'sw_down_fluxc_ext')
            except:
                if comm.rank == 0:
                    print('sw_down_fluxc_ext error')    
            try:
                io.output_profile(self.lw_tt_ext,'lw_tt_ext')
            except:
                if comm.rank == 0:
                    print('lw_tt_ext error')   
            try:
                io.output_profile(self.sw_tt_ext,'sw_tt_ext')
            except:
                if comm.rank == 0:
                    print('sw_tt_ext error')   
            try:
                io.output_profile(self.total_tt_ext,'total_tt_ext')
            except:
                if comm.rank == 0:
                    print('total_tt_ext error')  
            try:
                io.output_profile(self.lw_ttc_ext,'lw_ttc_ext')
            except:
                if comm.rank == 0:
                    print('lw_ttc_ext error')   
            try:
                io.output_profile(self.sw_ttc_ext,'sw_ttc_ext')
            except:
                if comm.rank == 0:
                    print('sw_ttc_ext error')   
            try:
                io.output_profile(self.total_ttc_ext,'total_ttc_ext')
            except:
                if comm.rank == 0:
                    print('total_ttc_ext error')  
            try:
                io.output_timeseries(self.tsfc_in_ext,'tsfc_in_ext')
            except:
                if comm.rank == 0:
                    print('tsfc_in_ext error')  
             
        return
