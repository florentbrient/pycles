#!python
#cython: boundscheck=False
#cython: wraparound=False
#cython: initializedcheck=False
#cython: cdivision=True

cimport Grid
cimport ReferenceState
cimport PrognosticVariables
cimport DiagnosticVariables
from NetCDFIO cimport NetCDFIO_Stats
cimport ParallelMPI


import numpy as np
cimport numpy as np
import sys
import netCDF4 as nc
from libc.math cimport pow, cbrt, exp
from thermodynamic_functions import exner, cpm
include 'parameters.pxi'

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


cdef class Radiation_RRTM:
    def __init__(self, namelist, ParallelMPI.ParallelMPI Pa):
        casename = namelist['radiation']['name']
        if casename == 'RRTM':
            self.scheme = RadiationRRTM()
        return
        
    cpdef initialize(self,case_dict,Grid.Grid Gr, ReferenceState.ReferenceState Ref,
                     NetCDFIO_Stats NS, ParallelMPI.ParallelMPI Pa):
        self.scheme.initialize(case_dict, Gr, Ref, NS, Pa)
        return

    cpdef update(self, Grid.Grid Gr, ReferenceState.ReferenceState Ref,
                 PrognosticVariables.PrognosticVariables PV, DiagnosticVariables.DiagnosticVariables DV,
                 ParallelMPI.ParallelMPI Pa):
        self.scheme.update(Gr, Ref, PV, DV, Pa)
        return

    cpdef stats_io(self, Grid.Grid Gr, ReferenceState.ReferenceState Ref,
                   PrognosticVariables.PrognosticVariables PV, DiagnosticVariables.DiagnosticVariables DV,
                   NetCDFIO_Stats NS, ParallelMPI.ParallelMPI Pa):
        self.scheme.stats_io(Gr, PV, DV, NS, Pa)
        return
        
cdef class RadiationRRTM:
    def __init__(self,Grid.Grid Gr):

        self.lw_flux_down = np.zeros(Gr.dims.nlg[2],dtype=np.double,order='c')
        self.lw_flux_up = np.zeros(Gr.dims.nlg[2],dtype=np.double,order='c')
        self.sw_flux_down = np.zeros(Gr.dims.nlg[2],dtype=np.double,order='c')
        self.sw_flux_up = np.zeros(Gr.dims.nlg[2],dtype=np.double,order='c')
        
        self.lw_flux_down_clear = np.zeros(Gr.dims.nlg[2],dtype=np.double,order='c')
        self.lw_flux_up_clear = np.zeros(Gr.dims.nlg[2],dtype=np.double,order='c')
        self.sw_flux_down_clear = np.zeros(Gr.dims.nlg[2],dtype=np.double,order='c')
        self.sw_flux_up_clear = np.zeros(Gr.dims.nlg[2],dtype=np.double,order='c')
        self.dTdt_clear = np.zeros(Gr.dims.nlg[2],dtype=np.double,order='c')
        
        self.dTdt_rad = np.zeros(Gr.dims.nlg[0]* Gr.dims.nlg[1]* Gr.dims.nlg[2], dtype=np.double, order='c')
        self.dsdt_rad = np.zeros(Gr.dims.nlg[0]* Gr.dims.nlg[1]* Gr.dims.nlg[2], dtype=np.double, order='c')
        
        self.osr = np.zeros(Gr.dims.nlg[0]* Gr.dims.nlg[1], dtype=np.double, order='c')
        self.olr = np.zeros(Gr.dims.nlg[0]* Gr.dims.nlg[1], dtype=np.double, order='c')
        
        self.rel_tropo_T = None
        self.rel_toa_T   = None
        self.n_ext = None
        
        # CGILS-specific: 
        self.ref_aloft = None
        #self.pressure_aloft = None
        #self.p0i_aloft = None
        #self.temperature_aloft = None
        #self.y_vapor_aloft = None
        #self.y_cond_aloft = None
        #self.cloud_fraction_aloft = None
        
        # Test: 
        #self.p0i_ext = None
        #self.pressure_ext = None
        #self.temperature_ext = None
        #self.y_vapor_ext = None
        #self.y_cond_ext = None
        #self.cloud_fraction_ext = None
        
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
        
        return
        
    cpdef initialize(self,case_dict,Grid.Grid Gr, ReferenceState.ReferenceState Ref,
                     NetCDFIO_Stats NS, ParallelMPI.ParallelMPI Pa):
        cdef :
             Py_ssize_t gw = Gr.dims.gw
             Py_ssize_t k
        	 
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
            #self.p0i_aloft[0] = basicstate.p0_i_global_ghosted[grid.gw + grid.nz-1]
            #self.p0i_aloft[0] = Ref.p0_i_global_ghosted[Grid.gw + Gr.dims.n[2]-1] # or ? extract_local_ghosted (fbrient)
            self.p0i_aloft[0] = Ref.p[Grid.gw + Gr.dims.n[2]-1] # or ? extract_local_ghosted (fbrient)

            for k in xrange(1,self.n_ext+1):
                self.p0i_aloft[k] = self.p0i_aloft[k-1] - self.p0i_aloft[0]/self.n_ext
                self.pressure_aloft[k-1] = (self.p0i_aloft[k] + self.p0i_aloft[k-1])*0.5

            p = case_dict['forcing']['pressure']
            # Interpolate for temperature and moisture aloft (above of domain top) 
            try:
                self.temperature_aloft = np.interp(self.pressure_aloft,p,case_dict['forcing']['t_ref'])
            except:
                self.temperature_aloft = np.interp(self.pressure_aloft,p,case_dict['forcing']['thl_ref']*exner(Ref.p0))#(p/basicstate.p00)**(Ra/cpa))
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
        pressure    = Ref.p0
        p0i         = Ref.p0_half
        self.pressure_ext  = np.zeros(Gr.dims.n[2] + self.n_ext   ,dtype=np.double,order='F')
        self.p0i_ext       = np.zeros(Gr.dims.n[2] + self.n_ext+1 ,dtype=np.double,order='F')       
        plev        = self.p0i_ext

        for k in xrange(Gr.dims.n[2] + self.n_ext):
            if k <= Gr.dims.n[2] :
                self.p0i_ext[k] = p0i[gw + k-1]
            elif self.ref_aloft:
                self.p0i_ext[k] = self.p0i_aloft[k-Gr.dims.n[2]]
            else:
                self.p0i_ext[k] = self.p0i_ext[k-1] - self.p0i_ext[Gr.dims.n[2]]/self.n_ext
        for k in xrange(Gr.dims.n[2] + self.n_ext):
            if k <= (Gr.dims.n[2]-1):
                self.pressure_ext[k] = pressure[k]
            else:
                self.pressure_ext[k] = (self.p0i_ext[k] + self.p0i_ext[k+1])*0.5     

        # Added: Initialize rrtmg_lw and rrtmg_sw
        cdef double cpdair = np.float64(cpm)
        c_rrtmg_lw_init(&cpdair)
        c_rrtmg_sw_init(&cpdair)
        
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
            
        lw_gas = nc.Dataset(lw_input_file, 'r+', format='NETCDF4')
        
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
        trpath = np.zeros((Gr.dims.n[2] + self.n_ext+1,9),dtype=np.double,order='F')
        for k in xrange(Gr.dims.n[2] + self.n_ext):
            plev[k]   = self.p0i_ext[k]/100.0
            
        for i in xrange(1,Gr.dims.n[2] + self.n_ext+1):
            trpath[i,:] = trpath[i-1,:]
            if (plev[i-1] > lw_pressure[0]):
                trpath[i,:] = trpath[i,:] + (plev[i-1] - np.max((plev[i],lw_pressure[0])))/g*trace[:,0]
            for m in xrange(1,self.lw_np):
                #print i, m
                plow = np.min((plev[i-1],np.max((plev[i], lw_pressure[m-1]))))
                pupp = np.min((plev[i-1],np.max((plev[i], lw_pressure[m]))))
                if (plow > pupp):
                    pmid   = 0.5*(plow+pupp)
                    wgtlow = (pmid-lw_pressure[m])/(lw_pressure[m-1]-lw_pressure[m])
                    wgtupp = (lw_pressure[m-1]-pmid)/(lw_pressure[m-1]-lw_pressure[m])
                    trpath[i,:] = trpath[i,:] + (plow-pupp)/g*(wgtlow*trace[:,m-1]  + wgtupp*trace[:,m])
            if (plev[i] < lw_pressure[self.lw_np-1]):
                trpath[i,:] = trpath[i,:] + (np.min((plev[i-1],lw_pressure[self.lw_np-1]))-plev[i])/g*trace[:,self.lw_np-1]

        tmpTrace = np.zeros((Gr.dims.n[2] + self.n_ext,9),dtype=np.double,order='F')
        
        self.o3vmr    = np.zeros(Gr.dims.n[2] + self.n_ext,dtype=np.double,order='F')
        self.co2vmr   = np.zeros(Gr.dims.n[2] + self.n_ext,dtype=np.double,order='F')
        self.ch4vmr   = np.zeros(Gr.dims.n[2] + self.n_ext,dtype=np.double,order='F')
        self.n2ovmr   = np.zeros(Gr.dims.n[2] + self.n_ext,dtype=np.double,order='F')
        self.o2vmr    = np.zeros(Gr.dims.n[2] + self.n_ext,dtype=np.double,order='F')
        self.cfc11vmr = np.zeros(Gr.dims.n[2] + self.n_ext,dtype=np.double,order='F')
        self.cfc12vmr = np.zeros(Gr.dims.n[2] + self.n_ext,dtype=np.double,order='F')
        self.cfc22vmr = np.zeros(Gr.dims.n[2] + self.n_ext,dtype=np.double,order='F')
        self.ccl4vmr  = np.zeros(Gr.dims.n[2] + self.n_ext,dtype=np.double,order='F')
        
        if self.use_o3in == False:
            self.o3vmr[:]  = tmpTrace[:,0]
        else:
            # o3_trace, o3_pressure
            trpath_o3 = np.zeros(Gr.dims.n[2] + self.n_ext+1,dtype=np.double,order='F')
            for k in xrange(Gr.dims.n[2] + self.n_ext):
                plev[k]   = self.p0i_ext[k]/100.0
                
            self.o3_np = o3_trace.shape[0]
            for i in xrange(1,Gr.dims.n[2] + self.n_ext+1):
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
            tmpTrace_o3 = np.zeros(Gr.dims.n[2] + self.n_ext,dtype=np.double,order='F')
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
#            type(thermodynamics.y_vapor)
            print('RadiationFMS: Thermodynamics With Vapor')
        except:
            self.is_vapor = False
            print('RadiationFMS: Thermodynamics WithOUT Vapor')
            
        self.is_liquid  = True
        try:
#            type(thermodynamics.y_liquid)
            print('RadiationFMS: Thermodynamics With Liquid')
        except:
            self.is_liquid = False   
            print('RadiationFMS: Thermodynamics WithOUT Liquid')
                     
        self.is_ice  = True
        try:
#            type(thermodynamics.y_ice)
            print('RadiationFMS: Thermodynamics With Ice')
        except:
            self.is_ice = False
            print('RadiationFMS: Thermodynamics WithOUT Ice')
            
        #Initialize statistical output file
#        self.init_output(Grid, io)
        
    cpdef stats_io(self,Grid.Grid Gr, ReferenceState.ReferenceState Ref,
                     NetCDFIO_Stats NS, ParallelMPI.ParallelMPI Pa):
            
        return            
            
        



        
        
        