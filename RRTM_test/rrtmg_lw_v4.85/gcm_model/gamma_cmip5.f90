      program gamma_cmip5
      implicit none
      
      ! Routine to calculate offline radiative fluxes from CMIP5
      ! You give cloud fraction profile from monthly data and calculate CRE ... (McCoy et. al?)
	  !       
      
      
      ! ----- Input -----
! Note: All volume mixing ratios are in dimensionless units of mole fraction obtained
! by scaling mass mixing ratio (g/g) with the appropriate molecular weights (g/mol) 
      integer(kind=im), intent(in) :: ncol            ! Number of horizontal columns
      integer(kind=im), intent(in) :: nlay            ! Number of model layers
      integer(kind=im), intent(inout) :: icld         ! Cloud overlap method
                                                      !    0: Clear only
                                                      !    1: Random
                                                      !    2: Maximum/random
                                                      !    3: Maximum
      integer(kind=im), intent(in) :: idrv            ! Flag for calculation of dFdT, the change
                                                      !    in upward flux as a function of 
                                                      !    surface temperature [0=off, 1=on]
                                                      !    0: Normal forward calculation
                                                      !    1: Normal forward calculation with
                                                      !       duflx_dt and duflxc_dt output

      real(kind=rb), intent(in) :: play(:,:)          ! Layer pressures (hPa, mb)
                                                      !    Dimensions: (ncol,nlay)
      real(kind=rb), intent(in) :: plev(:,:)          ! Interface pressures (hPa, mb)
                                                      !    Dimensions: (ncol,nlay+1)
      real(kind=rb), intent(in) :: tlay(:,:)          ! Layer temperatures (K)
                                                      !    Dimensions: (ncol,nlay)
      real(kind=rb), intent(in) :: tlev(:,:)          ! Interface temperatures (K)
                                                      !    Dimensions: (ncol,nlay+1)
      real(kind=rb), intent(in) :: tsfc(:)            ! Surface temperature (K)
                                                      !    Dimensions: (ncol)
      real(kind=rb), intent(in) :: h2ovmr(:,:)        ! H2O volume mixing ratio
                                                      !    Dimensions: (ncol,nlay)
      real(kind=rb), intent(in) :: o3vmr(:,:)         ! O3 volume mixing ratio
                                                      !    Dimensions: (ncol,nlay)
      real(kind=rb), intent(in) :: co2vmr(:,:)        ! CO2 volume mixing ratio
                                                      !    Dimensions: (ncol,nlay)
      real(kind=rb), intent(in) :: ch4vmr(:,:)        ! Methane volume mixing ratio
                                                      !    Dimensions: (ncol,nlay)
      real(kind=rb), intent(in) :: n2ovmr(:,:)        ! Nitrous oxide volume mixing ratio
                                                      !    Dimensions: (ncol,nlay)
      real(kind=rb), intent(in) :: o2vmr(:,:)         ! Oxygen volume mixing ratio
                                                      !    Dimensions: (ncol,nlay)
      real(kind=rb), intent(in) :: cfc11vmr(:,:)      ! CFC11 volume mixing ratio
                                                      !    Dimensions: (ncol,nlay)
      real(kind=rb), intent(in) :: cfc12vmr(:,:)      ! CFC12 volume mixing ratio
                                                      !    Dimensions: (ncol,nlay)
      real(kind=rb), intent(in) :: cfc22vmr(:,:)      ! CFC22 volume mixing ratio
                                                      !    Dimensions: (ncol,nlay)
      real(kind=rb), intent(in) :: ccl4vmr(:,:)       ! CCL4 volume mixing ratio
                                                      !    Dimensions: (ncol,nlay)
      real(kind=rb), intent(in) :: emis(:,:)          ! Surface emissivity
                                                      !    Dimensions: (ncol,nbndlw)

      integer(kind=im), intent(in) :: inflglw         ! Flag for cloud optical properties
      integer(kind=im), intent(in) :: iceflglw        ! Flag for ice particle specification
      integer(kind=im), intent(in) :: liqflglw        ! Flag for liquid droplet specification

      real(kind=rb), intent(in) :: cldfmcl(:,:,:)     ! Cloud fraction
                                                      !    Dimensions: (ngptlw,ncol,nlay)
      real(kind=rb), intent(in) :: ciwpmcl(:,:,:)     ! In-cloud ice water path (g/m2)
                                                      !    Dimensions: (ngptlw,ncol,nlay)
      real(kind=rb), intent(in) :: clwpmcl(:,:,:)     ! In-cloud liquid water path (g/m2)
                                                      !    Dimensions: (ngptlw,ncol,nlay)
      real(kind=rb), intent(in) :: reicmcl(:,:)       ! Cloud ice particle effective size (microns)
                                                      !    Dimensions: (ncol,nlay)
                                                      ! specific definition of reicmcl depends on setting of iceflglw:
                                                      ! iceflglw = 0: ice effective radius, r_ec, (Ebert and Curry, 1992),
                                                      !               r_ec must be >= 10.0 microns
                                                      ! iceflglw = 1: ice effective radius, r_ec, (Ebert and Curry, 1992),
                                                      !               r_ec range is limited to 13.0 to 130.0 microns
                                                      ! iceflglw = 2: ice effective radius, r_k, (Key, Streamer Ref. Manual, 1996)
                                                      !               r_k range is limited to 5.0 to 131.0 microns
                                                      ! iceflglw = 3: generalized effective size, dge, (Fu, 1996),
                                                      !               dge range is limited to 5.0 to 140.0 microns
                                                      !               [dge = 1.0315 * r_ec]
      real(kind=rb), intent(in) :: relqmcl(:,:)       ! Cloud water drop effective radius (microns)
                                                      !    Dimensions: (ncol,nlay)
      real(kind=rb), intent(in) :: taucmcl(:,:,:)     ! In-cloud optical depth
                                                      !    Dimensions: (ngptlw,ncol,nlay)
!      real(kind=rb), intent(in) :: ssacmcl(:,:,:)    ! In-cloud single scattering albedo
                                                      !    Dimensions: (ngptlw,ncol,nlay)
                                                      !   for future expansion
                                                      !   lw scattering not yet available
!      real(kind=rb), intent(in) :: asmcmcl(:,:,:)    ! In-cloud asymmetry parameter
                                                      !    Dimensions: (ngptlw,ncol,nlay)
                                                      !   for future expansion
                                                      !   lw scattering not yet available
      real(kind=rb), intent(in) :: tauaer(:,:,:)      ! aerosol optical depth
                                                      !   at mid-point of LW spectral bands
                                                      !    Dimensions: (ncol,nlay,nbndlw)
!      real(kind=rb), intent(in) :: ssaaer(:,:,:)     ! aerosol single scattering albedo
                                                      !    Dimensions: (ncol,nlay,nbndlw)
                                                      !   for future expansion 
                                                      !   (lw aerosols/scattering not yet available)
!      real(kind=rb), intent(in) :: asmaer(:,:,:)     ! aerosol asymmetry parameter
                                                      !    Dimensions: (ncol,nlay,nbndlw)
                                                      !   for future expansion 
                                                      !   (lw aerosols/scattering not yet available)

! ----- Output -----

      real(kind=rb), intent(out) :: uflx(:,:)         ! Total sky longwave upward flux (W/m2)
                                                      !    Dimensions: (ncol,nlay+1)
      real(kind=rb), intent(out) :: dflx(:,:)         ! Total sky longwave downward flux (W/m2)
                                                      !    Dimensions: (ncol,nlay+1)
      real(kind=rb), intent(out) :: hr(:,:)           ! Total sky longwave radiative heating rate (K/d)
                                                      !    Dimensions: (ncol,nlay)
      real(kind=rb), intent(out) :: uflxc(:,:)        ! Clear sky longwave upward flux (W/m2)
                                                      !    Dimensions: (ncol,nlay+1)
      real(kind=rb), intent(out) :: dflxc(:,:)        ! Clear sky longwave downward flux (W/m2)
                                                      !    Dimensions: (ncol,nlay+1)
      real(kind=rb), intent(out) :: hrc(:,:)          ! Clear sky longwave radiative heating rate (K/d)
                                                      !    Dimensions: (ncol,nlay)

! ----- Optional Output -----
      real(kind=rb), intent(out), optional :: duflx_dt(:,:)     
                                                      ! change in upward longwave flux (w/m2/K)
                                                      ! with respect to surface temperature
                                                      !    Dimensions: (ncol,nlay+1)
      real(kind=rb), intent(out), optional :: duflxc_dt(:,:)    
                                                      ! change in clear sky upward longwave flux (w/m2/K)
                                                      ! with respect to surface temperature
                                                      !    Dimensions: (ncol,nlay+1)
      
      ! from rrtmg_lw_rad.f90 (line 464)
	  call rrtmg_lw_ini(cpdair)
      
      
      
      ! call routine in rrtmg_lw_rad.f90
      call rrtmg_lw &
            (ncol    ,nlay    ,icld    ,idrv    , &
             play    ,plev    ,tlay    ,tlev    ,tsfc    , & 
             h2ovmr  ,o3vmr   ,co2vmr  ,ch4vmr  ,n2ovmr  ,o2vmr , &
             cfc11vmr,cfc12vmr,cfc22vmr,ccl4vmr ,emis    , &
             inflglw ,iceflglw,liqflglw,cldfmcl , &
             taucmcl ,ciwpmcl ,clwpmcl ,reicmcl ,relqmcl , &
             tauaer  , &
             uflx    ,dflx    ,hr      ,uflxc   ,dflxc,  hrc, &
             duflx_dt,duflxc_dt )
      
      
      
      
      
      
      end