cimport Grid
cimport ReferenceState
cimport PrognosticVariables
cimport DiagnosticVariables
from NetCDFIO cimport NetCDFIO_Stats
cimport ParallelMPI

cdef class Radiation_RRTM:
    cdef:
        object scheme
    cpdef initialize(self,case_dict,Grid.Grid Gr, ReferenceState.ReferenceState Ref,
                     NetCDFIO_Stats NS, ParallelMPI.ParallelMPI Pa)
    cpdef update(self, Grid.Grid Gr, ReferenceState.ReferenceState Ref,
                 PrognosticVariables.PrognosticVariables PV, DiagnosticVariables.DiagnosticVariables DV,
                 ParallelMPI.ParallelMPI Pa)
    cpdef stats_io(self, Grid.Grid Gr, ReferenceState.ReferenceState Ref,
                   PrognosticVariables.PrognosticVariables PV, DiagnosticVariables.DiagnosticVariables DV,
                   NetCDFIO_Stats NS, ParallelMPI.ParallelMPI Pa)

cdef class RadiationRRTM:
    cdef:
        double lw_flux_down
        double lw_flux_up
        double sw_flux_down
        double sw_flux_up
        
        double lw_flux_down_clear
        double lw_flux_up_clear
        double sw_flux_down_clear
        double sw_flux_up_clear
        double dTdt_clear
        
        double dTdt_rad
        double dsdt_rad
        
        double osr
        double olr
        
        double rel_tropo_T
        double rel_toa_T  
        double n_ext
        
        # CGILS-specific: 
        double ref_aloft
        double[:] pressure_aloft
        double[:] p0i_aloft
        double[:] temperature_aloft
        double[:] y_vapor_aloft
        double[:] y_cond_aloft
        double[:] cloud_fraction_aloft
        
        # Test: 
        double[:] p0i_ext
        double[:] pressure_ext
        double[:] temperature_ext
        double[:] y_vapor_ext
        double[:] y_cond_ext
        double[:] cloud_fraction_ext
        double[:] plev
        
        double is_vapor 
        double is_liquid
        double is_ice   
        
        # Added for PyRRTM
        double lw_input_file
        double lw_gas   
        double trace    
        double co2_factor
        double h2o_factor
        double lw_np
        double uniform_reliq
        double smooth_qc
        double radiation_frequency
        
        double o3vmr   
        double co2vmr  
        double ch4vmr  
        double n2ovmr  
        double o2vmr   
        double cfc11vmr
        double cfc12vmr
        double cfc22vmr
        double ccl4vmr 
        
        double dyofyr  
        double adjes   
        double scon    
        double coszen  
        double adif    
        
        # Test
        double play_in_ext 
        double plev_in_ext     
        double tlay_in_ext     
        double tlev_in_ext     
        double h2ovmr_in_ext 
        double cldfr_in_ext  
        double cliqwp_in_ext 
        double reliq_in_ext 
        double tsfc_in_ext  
         
        double debug_mode   
         
        

    cpdef initialize(self,case_dict,Grid.Grid Gr, ReferenceState.ReferenceState Ref,
                     NetCDFIO_Stats NS, ParallelMPI.ParallelMPI Pa)
                             
    cpdef stats_io(self,Grid.Grid Gr, ReferenceState.ReferenceState Ref,
                     NetCDFIO_Stats NS, ParallelMPI.ParallelMPI Pa)
                 
