!     path:      $Source: /storm/rc1/cvsroot/rc/rrtmg_lw/src/rrtmg_lw_rad.f90,v $
!     author:    $Author: miacono $
!     revision:  $Revision: 1.12 $
!     created:   $Date: 2011/04/08 20:25:01 $
!

       module rrtmg_lw_rad

!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2002-2009, Atmospheric & Environmental Research, Inc. (AER).  |
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
!  --------------------------------------------------------------------------
!
! ****************************************************************************
! *                                                                          *
! *                              RRTMG_LW                                    *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                   a rapid radiative transfer model                       *
! *                       for the longwave region                            * 
! *             for application to general circulation models                *
! *                                                                          *
! *                                                                          *
! *            Atmospheric and Environmental Research, Inc.                  *
! *                        131 Hartwell Avenue                               *
! *                        Lexington, MA 02421                               *
! *                                                                          *
! *                                                                          *
! *                           Eli J. Mlawer                                  *
! *                        Jennifer S. Delamere                              *
! *                         Michael J. Iacono                                *
! *                         Shepard A. Clough                                *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                       email:  miacono@aer.com                            *
! *                       email:  emlawer@aer.com                            *
! *                       email:  jdelamer@aer.com                           *
! *                                                                          *
! *        The authors wish to acknowledge the contributions of the          *
! *        following people:  Steven J. Taubman, Karen Cady-Pereira,         *
! *        Patrick D. Brown, Ronald E. Farren, Luke Chen, Robert Bergstrom.  *
! *                                                                          *
! ****************************************************************************

! -------- Modules --------
      use rrlw_vsn
      use mcica_subcol_gen_lw, only: mcica_subcol_lw
      use rrtmg_lw_cldprmc, only: cldprmc
! *** Move the required call to rrtmg_lw_ini below and the following 
! use association to the GCM initialization area ***
!      use rrtmg_lw_init, only: rrtmg_lw_ini
      use rrtmg_lw_rtrnmc, only: rtrnmc
      use rrtmg_lw_setcoef, only: setcoef
      use rrtmg_lw_taumol, only: taumol

      implicit none

! public interfaces/functions/subroutines
      public :: rrtmg_lw, inatm

!------------------------------------------------------------------
      contains
!------------------------------------------------------------------

!------------------------------------------------------------------
! Public subroutines
!------------------------------------------------------------------

      subroutine rrtmg_lw &
            (ncol    ,nlay    ,icld    ,idrv    , &
             play    ,plev    ,tlay    ,tlev    ,tsfc    , & 
             h2ovmr  ,o3vmr   ,co2vmr  ,ch4vmr  ,n2ovmr  ,o2vmr , &
             cfc11vmr,cfc12vmr,cfc22vmr,ccl4vmr ,emis    , &
             inflglw ,iceflglw,liqflglw,cldfmcl , &
             taucmcl ,ciwpmcl ,clwpmcl ,reicmcl ,relqmcl , &
             tauaer  , &
             uflx    ,dflx    ,hr      ,uflxc   ,dflxc,  hrc, &
             duflx_dt,duflxc_dt )

! -------- Description --------

! This program is the driver subroutine for RRTMG_LW, the AER LW radiation 
! model for application to GCMs, that has been adapted from RRTM_LW for
! improved efficiency.
!
! NOTE: The call to RRTMG_LW_INI should be moved to the GCM initialization
!  area, since this has to be called only once. 
!
! This routine:
!    a) calls INATM to read in the atmospheric profile from GCM;
!       all layering in RRTMG is ordered from surface to toa. 
!    b) calls CLDPRMC to set cloud optical depth for McICA based 
!       on input cloud properties 
!    c) calls SETCOEF to calculate various quantities needed for 
!       the radiative transfer algorithm
!    d) calls TAUMOL to calculate gaseous optical depths for each 
!       of the 16 spectral bands
!    e) calls RTRNMC (for both clear and cloudy profiles) to perform the
!       radiative transfer calculation using McICA, the Monte-Carlo 
!       Independent Column Approximation, to represent sub-grid scale 
!       cloud variability
!    f) passes the necessary fluxes and cooling rates back to GCM
!
! Two modes of operation are possible:
!     The mode is chosen by using either rrtmg_lw.nomcica.f90 (to not use
!     McICA) or rrtmg_lw.f90 (to use McICA) to interface with a GCM. 
!
!    1) Standard, single forward model calculation (imca = 0)
!    2) Monte Carlo Independent Column Approximation (McICA, Pincus et al., 
!       JC, 2003) method is applied to the forward model calculation (imca = 1)
!
! This call to RRTMG_LW must be preceeded by a call to the module
!     mcica_subcol_gen_lw.f90 to run the McICA sub-column cloud generator,
!     which will provide the cloud physical or cloud optical properties
!     on the RRTMG quadrature point (ngpt) dimension.
!     Two random number generators are available for use when imca = 1.
!     This is chosen by setting flag irnd on input to mcica_subcol_gen_lw.
!     1) KISSVEC (irnd = 0)
!     2) Mersenne-Twister (irnd = 1)
!
! Two methods of cloud property input are possible:
!     Cloud properties can be input in one of two ways (controlled by input 
!     flags inflglw, iceflglw, and liqflglw; see text file rrtmg_lw_instructions
!     and subroutine rrtmg_lw_cldprmc.f90 for further details):
!
!    1) Input cloud fraction and cloud optical depth directly (inflglw = 0)
!    2) Input cloud fraction and cloud physical properties (inflglw = 1 or 2);  
!       cloud optical properties are calculated by cldprmc or cldprmc based
!       on input settings of iceflglw and liqflglw.  Ice particle size provided
!       must be appropriately defined for the ice parameterization selected. 
!
! One method of aerosol property input is possible:
!     Aerosol properties can be input in only one way (controlled by input 
!     flag iaer; see text file rrtmg_lw_instructions for further details):
!
!    1) Input aerosol optical depth directly by layer and spectral band (iaer=10);
!       band average optical depth at the mid-point of each spectral band.
!       RRTMG_LW currently treats only aerosol absorption;
!       scattering capability is not presently available.
!
! The optional calculation of the change in upward flux as a function of surface 
! temperature is available (controlled by input flag idrv).  This can be utilized 
! to approximate adjustments to the upward flux profile caused only by a change in 
! surface temperature between full radiation calls.  This feature uses the pre-
! calculated derivative of the Planck function with respect to surface temperature. 
!
!    1) Normal forward calculation for the input profile (idrv=0)
!    2) Normal forward calculation with optional calculation of the change
!       in upward flux as a function of surface temperature for clear sky
!       and total sky flux.  Flux partial derivatives are provided in arrays
!       duflx_dt and duflxc_dt for total and clear sky.  (idrv=1)
!
!
! ------- Modifications -------
!
! This version of RRTMG_LW has been modified from RRTM_LW to use a reduced 
! set of g-points for application to GCMs.  
!
!-- Original version (derived from RRTM_LW), reduction of g-points, other
!   revisions for use with GCMs.  
!     1999: M. J. Iacono, AER, Inc.
!-- Adapted for use with NCAR/CAM.
!     May 2004: M. J. Iacono, AER, Inc.
!-- Revised to add McICA capability. 
!     Nov 2005: M. J. Iacono, AER, Inc.
!-- Conversion to F90 formatting for consistency with rrtmg_sw.
!     Feb 2007: M. J. Iacono, AER, Inc.
!-- Modifications to formatting to use assumed-shape arrays.
!     Aug 2007: M. J. Iacono, AER, Inc.
!-- Modified to add longwave aerosol absorption.
!     Apr 2008: M. J. Iacono, AER, Inc.
!-- Added capability to calculate derivative of upward flux wrt surface temperature. 
!     Nov 2009: M. J. Iacono, E. J. Mlawer, AER, Inc.

! --------- Modules ----------

      use parrrtm, only : nbndlw, ngptlw, maxxsec, mxmol
      use rrlw_con, only: fluxfac, heatfac, oneminus, pi
      use rrlw_wvn, only: ng, ngb, nspa, nspb, wavenum1, wavenum2, delwave

! ------- Declarations -------

! ----- Input -----
! Note: All volume mixing ratios are in dimensionless units of mole fraction obtained
! by scaling mass mixing ratio (g/g) with the appropriate molecular weights (g/mol) 
      integer(kind=4), intent(in) :: ncol            ! Number of horizontal columns
      integer(kind=4), intent(in) :: nlay            ! Number of model layers
      integer(kind=4), intent(inout) :: icld         ! Cloud overlap method
                                                      !    0: Clear only
                                                      !    1: Random
                                                      !    2: Maximum/random
                                                      !    3: Maximum
      integer(kind=4), intent(in) :: idrv            ! Flag for calculation of dFdT, the change
                                                      !    in upward flux as a function of 
                                                      !    surface temperature [0=off, 1=on]
                                                      !    0: Normal forward calculation
                                                      !    1: Normal forward calculation with
                                                      !       duflx_dt and duflxc_dt output

      real(kind=8), intent(in) :: play(:,:)          ! Layer pressures (hPa, mb)
                                                      !    Dimensions: (ncol,nlay)
      real(kind=8), intent(in) :: plev(:,:)          ! Interface pressures (hPa, mb)
                                                      !    Dimensions: (ncol,nlay+1)
      real(kind=8), intent(in) :: tlay(:,:)          ! Layer temperatures (K)
                                                      !    Dimensions: (ncol,nlay)
      real(kind=8), intent(in) :: tlev(:,:)          ! Interface temperatures (K)
                                                      !    Dimensions: (ncol,nlay+1)
      real(kind=8), intent(in) :: tsfc(:)            ! Surface temperature (K)
                                                      !    Dimensions: (ncol)
      real(kind=8), intent(in) :: h2ovmr(:,:)        ! H2O volume mixing ratio
                                                      !    Dimensions: (ncol,nlay)
      real(kind=8), intent(in) :: o3vmr(:,:)         ! O3 volume mixing ratio
                                                      !    Dimensions: (ncol,nlay)
      real(kind=8), intent(in) :: co2vmr(:,:)        ! CO2 volume mixing ratio
                                                      !    Dimensions: (ncol,nlay)
      real(kind=8), intent(in) :: ch4vmr(:,:)        ! Methane volume mixing ratio
                                                      !    Dimensions: (ncol,nlay)
      real(kind=8), intent(in) :: n2ovmr(:,:)        ! Nitrous oxide volume mixing ratio
                                                      !    Dimensions: (ncol,nlay)
      real(kind=8), intent(in) :: o2vmr(:,:)         ! Oxygen volume mixing ratio
                                                      !    Dimensions: (ncol,nlay)
      real(kind=8), intent(in) :: cfc11vmr(:,:)      ! CFC11 volume mixing ratio
                                                      !    Dimensions: (ncol,nlay)
      real(kind=8), intent(in) :: cfc12vmr(:,:)      ! CFC12 volume mixing ratio
                                                      !    Dimensions: (ncol,nlay)
      real(kind=8), intent(in) :: cfc22vmr(:,:)      ! CFC22 volume mixing ratio
                                                      !    Dimensions: (ncol,nlay)
      real(kind=8), intent(in) :: ccl4vmr(:,:)       ! CCL4 volume mixing ratio
                                                      !    Dimensions: (ncol,nlay)
      real(kind=8), intent(in) :: emis(:,:)          ! Surface emissivity
                                                      !    Dimensions: (ncol,nbndlw)

      integer(kind=4), intent(in) :: inflglw         ! Flag for cloud optical properties
      integer(kind=4), intent(in) :: iceflglw        ! Flag for ice particle specification
      integer(kind=4), intent(in) :: liqflglw        ! Flag for liquid droplet specification

      real(kind=8), intent(in) :: cldfmcl(:,:,:)     ! Cloud fraction
                                                      !    Dimensions: (ngptlw,ncol,nlay)
      real(kind=8), intent(in) :: ciwpmcl(:,:,:)     ! In-cloud ice water path (g/m2)
                                                      !    Dimensions: (ngptlw,ncol,nlay)
      real(kind=8), intent(in) :: clwpmcl(:,:,:)     ! In-cloud liquid water path (g/m2)
                                                      !    Dimensions: (ngptlw,ncol,nlay)
      real(kind=8), intent(in) :: reicmcl(:,:)       ! Cloud ice particle effective size (microns)
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
      real(kind=8), intent(in) :: relqmcl(:,:)       ! Cloud water drop effective radius (microns)
                                                      !    Dimensions: (ncol,nlay)
      real(kind=8), intent(in) :: taucmcl(:,:,:)     ! In-cloud optical depth
                                                      !    Dimensions: (ngptlw,ncol,nlay)
!      real(kind=8), intent(in) :: ssacmcl(:,:,:)    ! In-cloud single scattering albedo
                                                      !    Dimensions: (ngptlw,ncol,nlay)
                                                      !   for future expansion
                                                      !   lw scattering not yet available
!      real(kind=8), intent(in) :: asmcmcl(:,:,:)    ! In-cloud asymmetry parameter
                                                      !    Dimensions: (ngptlw,ncol,nlay)
                                                      !   for future expansion
                                                      !   lw scattering not yet available
      real(kind=8), intent(in) :: tauaer(:,:,:)      ! aerosol optical depth
                                                      !   at mid-point of LW spectral bands
                                                      !    Dimensions: (ncol,nlay,nbndlw)
!      real(kind=8), intent(in) :: ssaaer(:,:,:)     ! aerosol single scattering albedo
                                                      !    Dimensions: (ncol,nlay,nbndlw)
                                                      !   for future expansion 
                                                      !   (lw aerosols/scattering not yet available)
!      real(kind=8), intent(in) :: asmaer(:,:,:)     ! aerosol asymmetry parameter
                                                      !    Dimensions: (ncol,nlay,nbndlw)
                                                      !   for future expansion 
                                                      !   (lw aerosols/scattering not yet available)

! ----- Output -----

      real(kind=8), intent(out) :: uflx(:,:)         ! Total sky longwave upward flux (W/m2)
                                                      !    Dimensions: (ncol,nlay+1)
      real(kind=8), intent(out) :: dflx(:,:)         ! Total sky longwave downward flux (W/m2)
                                                      !    Dimensions: (ncol,nlay+1)
      real(kind=8), intent(out) :: hr(:,:)           ! Total sky longwave radiative heating rate (K/d)
                                                      !    Dimensions: (ncol,nlay)
      real(kind=8), intent(out) :: uflxc(:,:)        ! Clear sky longwave upward flux (W/m2)
                                                      !    Dimensions: (ncol,nlay+1)
      real(kind=8), intent(out) :: dflxc(:,:)        ! Clear sky longwave downward flux (W/m2)
                                                      !    Dimensions: (ncol,nlay+1)
      real(kind=8), intent(out) :: hrc(:,:)          ! Clear sky longwave radiative heating rate (K/d)
                                                      !    Dimensions: (ncol,nlay)

! ----- Optional Output -----
      real(kind=8), intent(out), optional :: duflx_dt(:,:)     
                                                      ! change in upward longwave flux (w/m2/K)
                                                      ! with respect to surface temperature
                                                      !    Dimensions: (ncol,nlay+1)
      real(kind=8), intent(out), optional :: duflxc_dt(:,:)    
                                                      ! change in clear sky upward longwave flux (w/m2/K)
                                                      ! with respect to surface temperature
                                                      !    Dimensions: (ncol,nlay+1)

! ----- Local -----

! Control
      integer(kind=4) :: nlayers             ! total number of layers
      integer(kind=4) :: istart              ! beginning band of calculation
      integer(kind=4) :: iend                ! ending band of calculation
      integer(kind=4) :: iout                ! output option flag (inactive)
      integer(kind=4) :: iaer                ! aerosol option flag
      integer(kind=4) :: iplon               ! column loop index
      integer(kind=4) :: imca                ! flag for mcica [0=off, 1=on]
      integer(kind=4) :: ims                 ! value for changing mcica permute seed
      integer(kind=4) :: k                   ! layer loop index
      integer(kind=4) :: ig                  ! g-point loop index

! Atmosphere
      real(kind=8) :: pavel(nlay+1)          ! layer pressures (mb) 
      real(kind=8) :: tavel(nlay+1)          ! layer temperatures (K)
      real(kind=8) :: pz(0:nlay+1)           ! level (interface) pressures (hPa, mb)
      real(kind=8) :: tz(0:nlay+1)           ! level (interface) temperatures (K)
      real(kind=8) :: tbound                 ! surface temperature (K)
      real(kind=8) :: coldry(nlay+1)         ! dry air column density (mol/cm2)
      real(kind=8) :: wbrodl(nlay+1)         ! broadening gas column density (mol/cm2)
      real(kind=8) :: wkl(mxmol,nlay+1)      ! molecular amounts (mol/cm-2)
      real(kind=8) :: wx(maxxsec,nlay+1)     ! cross-section amounts (mol/cm-2)
      real(kind=8) :: pwvcm                  ! precipitable water vapor (cm)
      real(kind=8) :: semiss(nbndlw)         ! lw surface emissivity
      real(kind=8) :: fracs(nlay+1,ngptlw)   ! 
      real(kind=8) :: taug(nlay+1,ngptlw)    ! gaseous optical depths
      real(kind=8) :: taut(nlay+1,ngptlw)    ! gaseous + aerosol optical depths

      real(kind=8) :: taua(nlay+1,nbndlw)    ! aerosol optical depth
!      real(kind=8) :: ssaa(nlay+1,nbndlw)   ! aerosol single scattering albedo
                                              !   for future expansion 
                                              !   (lw aerosols/scattering not yet available)
!      real(kind=8) :: asma(nlay+1,nbndlw)   ! aerosol asymmetry parameter
                                              !   for future expansion 
                                              !   (lw aerosols/scattering not yet available)

! Atmosphere - setcoef
      integer(kind=4) :: laytrop             ! tropopause layer index
      integer(kind=4) :: jp(nlay+1)          ! lookup table index 
      integer(kind=4) :: jt(nlay+1)          ! lookup table index 
      integer(kind=4) :: jt1(nlay+1)         ! lookup table index 
      real(kind=8) :: planklay(nlay+1,nbndlw)! 
      real(kind=8) :: planklev(0:nlay+1,nbndlw)! 
      real(kind=8) :: plankbnd(nbndlw)       ! 
      real(kind=8) :: dplankbnd_dt(nbndlw)   ! 

      real(kind=8) :: colh2o(nlay+1)         ! column amount (h2o)
      real(kind=8) :: colco2(nlay+1)         ! column amount (co2)
      real(kind=8) :: colo3(nlay+1)          ! column amount (o3)
      real(kind=8) :: coln2o(nlay+1)         ! column amount (n2o)
      real(kind=8) :: colco(nlay+1)          ! column amount (co)
      real(kind=8) :: colch4(nlay+1)         ! column amount (ch4)
      real(kind=8) :: colo2(nlay+1)          ! column amount (o2)
      real(kind=8) :: colbrd(nlay+1)         ! column amount (broadening gases)

      integer(kind=4) :: indself(nlay+1)
      integer(kind=4) :: indfor(nlay+1)
      real(kind=8) :: selffac(nlay+1)
      real(kind=8) :: selffrac(nlay+1)
      real(kind=8) :: forfac(nlay+1)
      real(kind=8) :: forfrac(nlay+1)

      integer(kind=4) :: indminor(nlay+1)
      real(kind=8) :: minorfrac(nlay+1)
      real(kind=8) :: scaleminor(nlay+1)
      real(kind=8) :: scaleminorn2(nlay+1)

      real(kind=8) :: &                      !
                         fac00(nlay+1), fac01(nlay+1), &
                         fac10(nlay+1), fac11(nlay+1) 
      real(kind=8) :: &                      !
                         rat_h2oco2(nlay+1),rat_h2oco2_1(nlay+1), &
                         rat_h2oo3(nlay+1),rat_h2oo3_1(nlay+1), &
                         rat_h2on2o(nlay+1),rat_h2on2o_1(nlay+1), &
                         rat_h2och4(nlay+1),rat_h2och4_1(nlay+1), &
                         rat_n2oco2(nlay+1),rat_n2oco2_1(nlay+1), &
                         rat_o3co2(nlay+1),rat_o3co2_1(nlay+1)

! Atmosphere/clouds - cldprop
      integer(kind=4) :: ncbands             ! number of cloud spectral bands
      integer(kind=4) :: inflag              ! flag for cloud property method
      integer(kind=4) :: iceflag             ! flag for ice cloud properties
      integer(kind=4) :: liqflag             ! flag for liquid cloud properties

! Atmosphere/clouds - cldprmc [mcica]
      real(kind=8) :: cldfmc(ngptlw,nlay+1)  ! cloud fraction [mcica]
      real(kind=8) :: ciwpmc(ngptlw,nlay+1)  ! in-cloud ice water path [mcica]
      real(kind=8) :: clwpmc(ngptlw,nlay+1)  ! in-cloud liquid water path [mcica]
      real(kind=8) :: relqmc(nlay+1)         ! liquid particle effective radius (microns)
      real(kind=8) :: reicmc(nlay+1)         ! ice particle effective size (microns)
      real(kind=8) :: taucmc(ngptlw,nlay+1)  ! in-cloud optical depth [mcica]
!      real(kind=8) :: ssacmc(ngptlw,nlay+1) ! in-cloud single scattering albedo [mcica]
                                              !   for future expansion 
                                              !   (lw scattering not yet available)
!      real(kind=8) :: asmcmc(ngptlw,nlay+1) ! in-cloud asymmetry parameter [mcica]
                                              !   for future expansion 
                                              !   (lw scattering not yet available)

! Output
      real(kind=8) :: totuflux(0:nlay+1)     ! upward longwave flux (w/m2)
      real(kind=8) :: totdflux(0:nlay+1)     ! downward longwave flux (w/m2)
      real(kind=8) :: fnet(0:nlay+1)         ! net longwave flux (w/m2)
      real(kind=8) :: htr(0:nlay+1)          ! longwave heating rate (k/day)
      real(kind=8) :: totuclfl(0:nlay+1)     ! clear sky upward longwave flux (w/m2)
      real(kind=8) :: totdclfl(0:nlay+1)     ! clear sky downward longwave flux (w/m2)
      real(kind=8) :: fnetc(0:nlay+1)        ! clear sky net longwave flux (w/m2)
      real(kind=8) :: htrc(0:nlay+1)         ! clear sky longwave heating rate (k/day)
      real(kind=8) :: dtotuflux_dt(0:nlay+1) ! change in upward longwave flux (w/m2/k)
                                              ! with respect to surface temperature
      real(kind=8) :: dtotuclfl_dt(0:nlay+1) ! change in clear sky upward longwave flux (w/m2/k)
                                              ! with respect to surface temperature

!
! Initializations

      oneminus = 1._8 - 1.e-6_8
      pi = 2._8 * asin(1._8)
      fluxfac = pi * 2.e4_8                  ! orig:   fluxfac = pi * 2.d4  
      istart = 1
      iend = 16
      iout = 0
      ims = 1

! Set imca to select calculation type:
!  imca = 0, use standard forward model calculation
!  imca = 1, use McICA for Monte Carlo treatment of sub-grid cloud variability

! *** This version uses McICA (imca = 1) ***

! Set icld to select of clear or cloud calculation and cloud overlap method  
! icld = 0, clear only
! icld = 1, with clouds using random cloud overlap
! icld = 2, with clouds using maximum/random cloud overlap
! icld = 3, with clouds using maximum cloud overlap (McICA only)
      if (icld.lt.0.or.icld.gt.3) icld = 2

! Set iaer to select aerosol option
! iaer = 0, no aerosols
! icld = 10, input total aerosol optical depth (tauaer) directly
      iaer = 10

! Call model and data initialization, compute lookup tables, perform
! reduction of g-points from 256 to 140 for input absorption coefficient 
! data and other arrays.
!
! In a GCM this call should be placed in the model initialization
! area, since this has to be called only once.  
!      call rrtmg_lw_ini(cpdair)

!  This is the main longitude/column loop within RRTMG.
      do iplon = 1, ncol

!  Prepare atmospheric profile from GCM for use in RRTMG, and define
!  other input parameters.  

         call inatm (iplon, nlay, icld, iaer, &
              play, plev, tlay, tlev, tsfc, h2ovmr, &
              o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr, cfc11vmr, cfc12vmr, &
              cfc22vmr, ccl4vmr, emis, inflglw, iceflglw, liqflglw, &
              cldfmcl, taucmcl, ciwpmcl, clwpmcl, reicmcl, relqmcl, tauaer, &
              nlayers, pavel, pz, tavel, tz, tbound, semiss, coldry, &
              wkl, wbrodl, wx, pwvcm, inflag, iceflag, liqflag, &
              cldfmc, taucmc, ciwpmc, clwpmc, reicmc, relqmc, taua)

!  For cloudy atmosphere, use cldprmc to set cloud optical properties based on
!  input cloud physical properties.  Select method based on choices described
!  in cldprmc.  Cloud fraction, water path, liquid droplet and ice particle
!  effective radius must be passed into cldprmc.  Cloud fraction and cloud
!  optical depth are transferred to rrtmg_lw arrays in cldprmc.  

         call cldprmc(nlayers, inflag, iceflag, liqflag, cldfmc, ciwpmc, &
                      clwpmc, reicmc, relqmc, ncbands, taucmc)

! Calculate information needed by the radiative transfer routine
! that is specific to this atmosphere, especially some of the 
! coefficients and indices needed to compute the optical depths
! by interpolating data from stored reference atmospheres. 

         call setcoef(nlayers, istart, pavel, tavel, tz, tbound, semiss, &
                      coldry, wkl, wbrodl, &
                      laytrop, jp, jt, jt1, planklay, planklev, plankbnd, &
                      idrv, dplankbnd_dt, &
                      colh2o, colco2, colo3, coln2o, colco, colch4, colo2, &
                      colbrd, fac00, fac01, fac10, fac11, &
                      rat_h2oco2, rat_h2oco2_1, rat_h2oo3, rat_h2oo3_1, &
                      rat_h2on2o, rat_h2on2o_1, rat_h2och4, rat_h2och4_1, &
                      rat_n2oco2, rat_n2oco2_1, rat_o3co2, rat_o3co2_1, &
                      selffac, selffrac, indself, forfac, forfrac, indfor, &
                      minorfrac, scaleminor, scaleminorn2, indminor)

!  Calculate the gaseous optical depths and Planck fractions for 
!  each longwave spectral band.

         call taumol(nlayers, pavel, wx, coldry, &
                     laytrop, jp, jt, jt1, planklay, planklev, plankbnd, &
                     colh2o, colco2, colo3, coln2o, colco, colch4, colo2, &
                     colbrd, fac00, fac01, fac10, fac11, &
                     rat_h2oco2, rat_h2oco2_1, rat_h2oo3, rat_h2oo3_1, &
                     rat_h2on2o, rat_h2on2o_1, rat_h2och4, rat_h2och4_1, &
                     rat_n2oco2, rat_n2oco2_1, rat_o3co2, rat_o3co2_1, &
                     selffac, selffrac, indself, forfac, forfrac, indfor, &
                     minorfrac, scaleminor, scaleminorn2, indminor, &
                     fracs, taug)


! Combine gaseous and aerosol optical depths, if aerosol active
         if (iaer .eq. 0) then
            do k = 1, nlayers
               do ig = 1, ngptlw
                  taut(k,ig) = taug(k,ig)
               enddo
            enddo
         elseif (iaer .eq. 10) then
            do k = 1, nlayers
               do ig = 1, ngptlw
                  taut(k,ig) = taug(k,ig) + taua(k,ngb(ig))
               enddo
            enddo
         endif

! Call the radiative transfer routine.
! Either routine can be called to do clear sky calculation.  If clouds
! are present, then select routine based on cloud overlap assumption
! to be used.  Clear sky calculation is done simultaneously.
! For McICA, RTRNMC is called for clear and cloudy calculations.

         call rtrnmc(nlayers, istart, iend, iout, pz, semiss, ncbands, &
                     cldfmc, taucmc, planklay, planklev, plankbnd, &
                     pwvcm, fracs, taut, &
                     totuflux, totdflux, fnet, htr, &
                     totuclfl, totdclfl, fnetc, htrc, &
                     idrv, dplankbnd_dt, dtotuflux_dt, dtotuclfl_dt )

!  Transfer up and down fluxes and heating rate to output arrays.
!  Vertical indexing goes from bottom to top; reverse here for GCM if necessary.

         do k = 0, nlayers
            uflx(iplon,k+1) = totuflux(k)
            dflx(iplon,k+1) = totdflux(k)
            uflxc(iplon,k+1) = totuclfl(k)
            dflxc(iplon,k+1) = totdclfl(k)
         enddo
         do k = 0, nlayers-1
            hr(iplon,k+1) = htr(k)
            hrc(iplon,k+1) = htrc(k)
         enddo

!  If idrv=1 option is active, transfer upward flux derivatives to output arrays.
         
         if (idrv .eq. 1) then 
            do k = 0, nlayers
               duflx_dt(iplon,k+1) = dtotuflux_dt(k)
               duflxc_dt(iplon,k+1) = dtotuclfl_dt(k)
            enddo
         endif

! End longitude/column loop
      enddo

      end subroutine rrtmg_lw

!***************************************************************************
      subroutine inatm (iplon, nlay, icld, iaer, &
              play, plev, tlay, tlev, tsfc, h2ovmr, &
              o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr, cfc11vmr, cfc12vmr, &
              cfc22vmr, ccl4vmr, emis, inflglw, iceflglw, liqflglw, &
              cldfmcl, taucmcl, ciwpmcl, clwpmcl, reicmcl, relqmcl, tauaer, &
              nlayers, pavel, pz, tavel, tz, tbound, semiss, coldry, &
              wkl, wbrodl, wx, pwvcm, inflag, iceflag, liqflag, &
              cldfmc, taucmc, ciwpmc, clwpmc, reicmc, relqmc, taua)
!***************************************************************************
!
!  Input atmospheric profile from GCM, and prepare it for use in RRTMG_LW.
!  Set other RRTMG_LW input parameters.  
!
!***************************************************************************

! --------- Modules ----------

      use parrrtm, only : nbndlw, ngptlw, nmol, maxxsec, mxmol
      use rrlw_con, only: fluxfac, heatfac, oneminus, pi, grav, avogad
      use rrlw_wvn, only: ng, nspa, nspb, wavenum1, wavenum2, delwave, ixindx

! ------- Declarations -------

! ----- Input -----
! Note: All volume mixing ratios are in dimensionless units of mole fraction obtained
! by scaling mass mixing ratio (g/g) with the appropriate molecular weights (g/mol) 
      integer(kind=4), intent(in) :: iplon           ! column loop index
      integer(kind=4), intent(in) :: nlay            ! Number of model layers
      integer(kind=4), intent(in) :: icld            ! clear/cloud and cloud overlap flag
      integer(kind=4), intent(in) :: iaer            ! aerosol option flag

      real(kind=8), intent(in) :: play(:,:)          ! Layer pressures (hPa, mb)
                                                      !    Dimensions: (ncol,nlay)
      real(kind=8), intent(in) :: plev(:,:)          ! Interface pressures (hPa, mb)
                                                      !    Dimensions: (ncol,nlay+1)
      real(kind=8), intent(in) :: tlay(:,:)          ! Layer temperatures (K)
                                                      !    Dimensions: (ncol,nlay)
      real(kind=8), intent(in) :: tlev(:,:)          ! Interface temperatures (K)
                                                      !    Dimensions: (ncol,nlay+1)
      real(kind=8), intent(in) :: tsfc(:)            ! Surface temperature (K)
                                                      !    Dimensions: (ncol)
      real(kind=8), intent(in) :: h2ovmr(:,:)        ! H2O volume mixing ratio
                                                      !    Dimensions: (ncol,nlay)
      real(kind=8), intent(in) :: o3vmr(:,:)         ! O3 volume mixing ratio
                                                      !    Dimensions: (ncol,nlay)
      real(kind=8), intent(in) :: co2vmr(:,:)        ! CO2 volume mixing ratio
                                                      !    Dimensions: (ncol,nlay)
      real(kind=8), intent(in) :: ch4vmr(:,:)        ! Methane volume mixing ratio
                                                      !    Dimensions: (ncol,nlay)
      real(kind=8), intent(in) :: n2ovmr(:,:)        ! Nitrous oxide volume mixing ratio
                                                      !    Dimensions: (ncol,nlay)
      real(kind=8), intent(in) :: o2vmr(:,:)         ! Oxygen volume mixing ratio
                                                      !    Dimensions: (ncol,nlay)
      real(kind=8), intent(in) :: cfc11vmr(:,:)      ! CFC11 volume mixing ratio
                                                      !    Dimensions: (ncol,nlay)
      real(kind=8), intent(in) :: cfc12vmr(:,:)      ! CFC12 volume mixing ratio
                                                      !    Dimensions: (ncol,nlay)
      real(kind=8), intent(in) :: cfc22vmr(:,:)      ! CFC22 volume mixing ratio
                                                      !    Dimensions: (ncol,nlay)
      real(kind=8), intent(in) :: ccl4vmr(:,:)       ! CCL4 volume mixing ratio
                                                      !    Dimensions: (ncol,nlay)
      real(kind=8), intent(in) :: emis(:,:)          ! Surface emissivity
                                                      !    Dimensions: (ncol,nbndlw)

      integer(kind=4), intent(in) :: inflglw         ! Flag for cloud optical properties
      integer(kind=4), intent(in) :: iceflglw        ! Flag for ice particle specification
      integer(kind=4), intent(in) :: liqflglw        ! Flag for liquid droplet specification

      real(kind=8), intent(in) :: cldfmcl(:,:,:)     ! Cloud fraction
                                                      !    Dimensions: (ngptlw,ncol,nlay)
      real(kind=8), intent(in) :: ciwpmcl(:,:,:)     ! In-cloud ice water path (g/m2)
                                                      !    Dimensions: (ngptlw,ncol,nlay)
      real(kind=8), intent(in) :: clwpmcl(:,:,:)     ! In-cloud liquid water path (g/m2)
                                                      !    Dimensions: (ngptlw,ncol,nlay)
      real(kind=8), intent(in) :: relqmcl(:,:)       ! Cloud water drop effective radius (microns)
                                                      !    Dimensions: (ncol,nlay)
      real(kind=8), intent(in) :: reicmcl(:,:)       ! Cloud ice effective size (microns)
                                                      !    Dimensions: (ncol,nlay)
      real(kind=8), intent(in) :: taucmcl(:,:,:)     ! In-cloud optical depth
                                                      !    Dimensions: (ngptlw,ncol,nlay)
      real(kind=8), intent(in) :: tauaer(:,:,:)      ! Aerosol optical depth
                                                      !    Dimensions: (ncol,nlay,nbndlw)

! ----- Output -----
! Atmosphere
      integer(kind=4), intent(out) :: nlayers        ! number of layers

      real(kind=8), intent(out) :: pavel(:)          ! layer pressures (mb) 
                                                      !    Dimensions: (nlay)
      real(kind=8), intent(out) :: tavel(:)          ! layer temperatures (K)
                                                      !    Dimensions: (nlay)
      real(kind=8), intent(out) :: pz(0:)            ! level (interface) pressures (hPa, mb)
                                                      !    Dimensions: (0:nlay)
      real(kind=8), intent(out) :: tz(0:)            ! level (interface) temperatures (K)
                                                      !    Dimensions: (0:nlay)
      real(kind=8), intent(out) :: tbound            ! surface temperature (K)
      real(kind=8), intent(out) :: coldry(:)         ! dry air column density (mol/cm2)
                                                      !    Dimensions: (nlay)
      real(kind=8), intent(out) :: wbrodl(:)         ! broadening gas column density (mol/cm2)
                                                      !    Dimensions: (nlay)
      real(kind=8), intent(out) :: wkl(:,:)          ! molecular amounts (mol/cm-2)
                                                      !    Dimensions: (mxmol,nlay)
      real(kind=8), intent(out) :: wx(:,:)           ! cross-section amounts (mol/cm-2)
                                                      !    Dimensions: (maxxsec,nlay)
      real(kind=8), intent(out) :: pwvcm             ! precipitable water vapor (cm)
      real(kind=8), intent(out) :: semiss(:)         ! lw surface emissivity
                                                      !    Dimensions: (nbndlw)

! Atmosphere/clouds - cldprmc
      integer(kind=4), intent(out) :: inflag         ! flag for cloud property method
      integer(kind=4), intent(out) :: iceflag        ! flag for ice cloud properties
      integer(kind=4), intent(out) :: liqflag        ! flag for liquid cloud properties

      real(kind=8), intent(out) :: cldfmc(:,:)       ! cloud fraction [mcica]
                                                      !    Dimensions: (ngptlw,nlay)
      real(kind=8), intent(out) :: ciwpmc(:,:)       ! in-cloud ice water path [mcica]
                                                      !    Dimensions: (ngptlw,nlay)
      real(kind=8), intent(out) :: clwpmc(:,:)       ! in-cloud liquid water path [mcica]
                                                      !    Dimensions: (ngptlw,nlay)
      real(kind=8), intent(out) :: relqmc(:)         ! liquid particle effective radius (microns)
                                                      !    Dimensions: (nlay)
      real(kind=8), intent(out) :: reicmc(:)         ! ice particle effective size (microns)
                                                      !    Dimensions: (nlay)
      real(kind=8), intent(out) :: taucmc(:,:)       ! in-cloud optical depth [mcica]
                                                      !    Dimensions: (ngptlw,nlay)
      real(kind=8), intent(out) :: taua(:,:)         ! aerosol optical depth
                                                      !    Dimensions: (nlay,nbndlw)


! ----- Local -----
      real(kind=8), parameter :: amd = 28.9660_8    ! Effective molecular weight of dry air (g/mol)
      real(kind=8), parameter :: amw = 18.0160_8    ! Molecular weight of water vapor (g/mol)
!      real(kind=8), parameter :: amc = 44.0098_8    ! Molecular weight of carbon dioxide (g/mol)
!      real(kind=8), parameter :: amo = 47.9998_8    ! Molecular weight of ozone (g/mol)
!      real(kind=8), parameter :: amo2 = 31.9999_8   ! Molecular weight of oxygen (g/mol)
!      real(kind=8), parameter :: amch4 = 16.0430_8  ! Molecular weight of methane (g/mol)
!      real(kind=8), parameter :: amn2o = 44.0128_8  ! Molecular weight of nitrous oxide (g/mol)
!      real(kind=8), parameter :: amc11 = 137.3684_8 ! Molecular weight of CFC11 (g/mol) - CCL3F
!      real(kind=8), parameter :: amc12 = 120.9138_8 ! Molecular weight of CFC12 (g/mol) - CCL2F2
!      real(kind=8), parameter :: amc22 = 86.4688_8  ! Molecular weight of CFC22 (g/mol) - CHCLF2
!      real(kind=8), parameter :: amcl4 = 153.823_8  ! Molecular weight of CCL4 (g/mol) - CCL4

! Set molecular weight ratios (for converting mmr to vmr)
!  e.g. h2ovmr = h2ommr * amdw)
      real(kind=8), parameter :: amdw = 1.607793_8  ! Molecular weight of dry air / water vapor
      real(kind=8), parameter :: amdc = 0.658114_8  ! Molecular weight of dry air / carbon dioxide
      real(kind=8), parameter :: amdo = 0.603428_8  ! Molecular weight of dry air / ozone
      real(kind=8), parameter :: amdm = 1.805423_8  ! Molecular weight of dry air / methane
      real(kind=8), parameter :: amdn = 0.658090_8  ! Molecular weight of dry air / nitrous oxide
      real(kind=8), parameter :: amdo2 = 0.905140_8 ! Molecular weight of dry air / oxygen
      real(kind=8), parameter :: amdc1 = 0.210852_8 ! Molecular weight of dry air / CFC11
      real(kind=8), parameter :: amdc2 = 0.239546_8 ! Molecular weight of dry air / CFC12

      integer(kind=4) :: isp, l, ix, n, imol, ib, ig   ! Loop indices
      real(kind=8) :: amm, amttl, wvttl, wvsh, summol  

! Add one to nlayers here to include extra model layer at top of atmosphere
      nlayers = nlay

!  Initialize all molecular amounts and cloud properties to zero here, then pass input amounts
!  into RRTM arrays below.

      wkl(:,:) = 0.0_8
      wx(:,:) = 0.0_8
      cldfmc(:,:) = 0.0_8
      taucmc(:,:) = 0.0_8
      ciwpmc(:,:) = 0.0_8
      clwpmc(:,:) = 0.0_8
      reicmc(:) = 0.0_8
      relqmc(:) = 0.0_8
      taua(:,:) = 0.0_8
      amttl = 0.0_8
      wvttl = 0.0_8
 
!  Set surface temperature.
      tbound = tsfc(iplon)

!  Install input GCM arrays into RRTMG_LW arrays for pressure, temperature,
!  and molecular amounts.  
!  Pressures are input in mb, or are converted to mb here.
!  Molecular amounts are input in volume mixing ratio, or are converted from 
!  mass mixing ratio (or specific humidity for h2o) to volume mixing ratio
!  here. These are then converted to molecular amount (molec/cm2) below.  
!  The dry air column COLDRY (in molec/cm2) is calculated from the level 
!  pressures, pz (in mb), based on the hydrostatic equation and includes a 
!  correction to account for h2o in the layer.  The molecular weight of moist 
!  air (amm) is calculated for each layer.  
!  Note: In RRTMG, layer indexing goes from bottom to top, and coding below
!  assumes GCM input fields are also bottom to top. Input layer indexing
!  from GCM fields should be reversed here if necessary.

      pz(0) = plev(iplon,1)
      tz(0) = tlev(iplon,1)
      do l = 1, nlayers
         pavel(l) = play(iplon,l)
         tavel(l) = tlay(iplon,l)
         pz(l) = plev(iplon,l+1)
         tz(l) = tlev(iplon,l+1)
! For h2o input in vmr:
         wkl(1,l) = h2ovmr(iplon,l)
! For h2o input in mmr:
!         wkl(1,l) = h2o(iplon,l)*amdw
! For h2o input in specific humidity;
!         wkl(1,l) = (h2o(iplon,l)/(1._8 - h2o(iplon,l)))*amdw
         wkl(2,l) = co2vmr(iplon,l)
         wkl(3,l) = o3vmr(iplon,l)
         wkl(4,l) = n2ovmr(iplon,l)
         wkl(6,l) = ch4vmr(iplon,l)
         wkl(7,l) = o2vmr(iplon,l)
         amm = (1._8 - wkl(1,l)) * amd + wkl(1,l) * amw            
         coldry(l) = (pz(l-1)-pz(l)) * 1.e3_8 * avogad / &
                     (1.e2_8 * grav * amm * (1._8 + wkl(1,l)))
      enddo

! Set cross section molecule amounts from input; convert to vmr if necessary
      do l=1, nlayers
         wx(1,l) = ccl4vmr(iplon,l)
         wx(2,l) = cfc11vmr(iplon,l)
         wx(3,l) = cfc12vmr(iplon,l)
         wx(4,l) = cfc22vmr(iplon,l)
      enddo      

! The following section can be used to set values for an additional layer (from
! the GCM top level to 1.e-4 mb) for improved calculation of TOA fluxes. 
! Temperature and molecular amounts in the extra model layer are set to 
! their values in the top GCM model layer, though these can be modified
! here if necessary. 
! If this feature is utilized, increase nlayers by one above, limit the two
! loops above to (nlayers-1), and set the top most (extra) layer values here. 

!      pavel(nlayers) = 0.5_8 * pz(nlayers-1)
!      tavel(nlayers) = tavel(nlayers-1)
!      pz(nlayers) = 1.e-4_8
!      tz(nlayers-1) = 0.5_8 * (tavel(nlayers)+tavel(nlayers-1))
!      tz(nlayers) = tz(nlayers-1)
!      wkl(1,nlayers) = wkl(1,nlayers-1)
!      wkl(2,nlayers) = wkl(2,nlayers-1)
!      wkl(3,nlayers) = wkl(3,nlayers-1)
!      wkl(4,nlayers) = wkl(4,nlayers-1)
!      wkl(6,nlayers) = wkl(6,nlayers-1)
!      wkl(7,nlayers) = wkl(7,nlayers-1)
!      amm = (1._8 - wkl(1,nlayers-1)) * amd + wkl(1,nlayers-1) * amw
!      coldry(nlayers) = (pz(nlayers-1)) * 1.e3_8 * avogad / &
!                        (1.e2_8 * grav * amm * (1._8 + wkl(1,nlayers-1)))
!      wx(1,nlayers) = wx(1,nlayers-1)
!      wx(2,nlayers) = wx(2,nlayers-1)
!      wx(3,nlayers) = wx(3,nlayers-1)
!      wx(4,nlayers) = wx(4,nlayers-1)

! At this point all molecular amounts in wkl and wx are in volume mixing ratio; 
! convert to molec/cm2 based on coldry for use in rrtm.  also, compute precipitable
! water vapor for diffusivity angle adjustments in rtrn and rtrnmr.

      do l = 1, nlayers
         summol = 0.0_8
         do imol = 2, nmol
            summol = summol + wkl(imol,l)
         enddo
         wbrodl(l) = coldry(l) * (1._8 - summol)
         do imol = 1, nmol
            wkl(imol,l) = coldry(l) * wkl(imol,l)
         enddo
         amttl = amttl + coldry(l)+wkl(1,l)
         wvttl = wvttl + wkl(1,l)
         do ix = 1,maxxsec
            if (ixindx(ix) .ne. 0) then
               wx(ixindx(ix),l) = coldry(l) * wx(ix,l) * 1.e-20_8
            endif
         enddo
      enddo

      wvsh = (amw * wvttl) / (amd * amttl)
      pwvcm = wvsh * (1.e3_8 * pz(0)) / (1.e2_8 * grav)

! Set spectral surface emissivity for each longwave band.  

      do n=1,nbndlw
         semiss(n) = emis(iplon,n)
!          semiss(n) = 1.0_8
      enddo

! Transfer aerosol optical properties to RRTM variable;
! modify to reverse layer indexing here if necessary.

     if (iaer .ge. 1) then
        do l = 1, nlayers
           do ib = 1, nbndlw
              taua(l,ib) = tauaer(iplon,l,ib)
           enddo
        enddo
      endif

! Transfer cloud fraction and cloud optical properties to RRTM variables,
! modify to reverse layer indexing here if necessary.

      if (icld .ge. 1) then 
         inflag = inflglw
         iceflag = iceflglw
         liqflag = liqflglw

! Move incoming GCM cloud arrays to RRTMG cloud arrays.
! For GCM input, incoming reicmcl is defined based on selected ice parameterization (inflglw)

         do l = 1, nlayers
            do ig = 1, ngptlw
               cldfmc(ig,l) = cldfmcl(ig,iplon,l)
               taucmc(ig,l) = taucmcl(ig,iplon,l)
               ciwpmc(ig,l) = ciwpmcl(ig,iplon,l)
               clwpmc(ig,l) = clwpmcl(ig,iplon,l)
            enddo
            reicmc(l) = reicmcl(iplon,l)
            relqmc(l) = relqmcl(iplon,l)
         enddo

! If an extra layer is being used in RRTMG, set all cloud properties to zero in the extra layer.

!         cldfmc(:,nlayers) = 0.0_8
!         taucmc(:,nlayers) = 0.0_8
!         ciwpmc(:,nlayers) = 0.0_8
!         clwpmc(:,nlayers) = 0.0_8
!         reicmc(nlayers) = 0.0_8
!         relqmc(nlayers) = 0.0_8
!         taua(nlayers,:) = 0.0_8

      endif
      
      end subroutine inatm

      end module rrtmg_lw_rad
