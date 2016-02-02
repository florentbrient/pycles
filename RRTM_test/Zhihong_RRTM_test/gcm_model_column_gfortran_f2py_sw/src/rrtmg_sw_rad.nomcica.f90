!     path:      $Source$
!     author:    $Author: miacono $
!     revision:  $Revision: 23308 $
!     created:   $Date: 2013-12-27 17:23:51 -0500 (Fri, 27 Dec 2013) $
!

       module rrtmg_sw_rad

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
! *                             RRTMG_SW                                     *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                 a rapid radiative transfer model                         *
! *                  for the solar spectral region                           *
! *           for application to general circulation models                  *
! *                                                                          *
! *                                                                          *
! *           Atmospheric and Environmental Research, Inc.                   *
! *                       131 Hartwell Avenue                                *
! *                       Lexington, MA 02421                                *
! *                                                                          *
! *                                                                          *
! *                          Eli J. Mlawer                                   *
! *                       Jennifer S. Delamere                               *
! *                        Michael J. Iacono                                 *
! *                        Shepard A. Clough                                 *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                                                                          *
! *                      email:  emlawer@aer.com                             *
! *                      email:  jdelamer@aer.com                            *
! *                      email:  miacono@aer.com                             *
! *                                                                          *
! *       The authors wish to acknowledge the contributions of the           *
! *       following people:  Steven J. Taubman, Patrick D. Brown,            *
! *       Ronald E. Farren, Luke Chen, Robert Bergstrom.                     *
! *                                                                          *
! ****************************************************************************

! --------- Modules ---------
      
      use rrsw_vsn
      use rrtmg_sw_cldprop, only: cldprop_sw
! *** Move the required call to rrtmg_sw_ini below and the following 
! use association to GCM initialization area ***
!      use rrtmg_sw_init, only: rrtmg_sw_ini
      use rrtmg_sw_setcoef, only: setcoef_sw
      use rrtmg_sw_spcvrt, only: spcvrt_sw

      implicit none

! public interfaces/functions/subroutines
      public :: rrtmg_sw, inatm_sw, earth_sun

!------------------------------------------------------------------
      contains
!------------------------------------------------------------------

!------------------------------------------------------------------
! Public subroutines
!------------------------------------------------------------------

      subroutine rrtmg_sw &
            (ncol    ,nlay    ,icld    ,iaer    ,idelm, &    ! idelm added by ZTAN
             play    ,plev    ,tlay    ,tlev    ,tsfc    , &
             h2ovmr  ,o3vmr   ,co2vmr  ,ch4vmr  ,n2ovmr  ,o2vmr, &
             asdir   ,asdif   ,aldir   ,aldif   , &
             coszen  ,adjes   ,dyofyr  ,scon    , &
             inflgsw ,iceflgsw,liqflgsw,cldfr   , &
             taucld  ,ssacld  ,asmcld  ,fsfcld  , &
             cicewp  ,cliqwp  ,reice   ,reliq   , &
             tauaer  ,ssaaer  ,asmaer  ,ecaer   , &
             swuflx  ,swdflx  ,swhr    ,swuflxc ,swdflxc ,swhrc, &
             dirdflux, difdflux)

! ------- Description -------

! This program is the driver for RRTMG_SW, the AER SW radiation model for 
!  application to GCMs, that has been adapted from RRTM_SW for improved
!  efficiency and to provide fractional cloudiness and cloud overlap
!  capability using McICA.
!
! Note: The call to RRTMG_SW_INI should be moved to the GCM initialization 
!  area, since this has to be called only once. 
!
! This routine
!    b) calls INATM_SW to read in the atmospheric profile from GCM;
!       all layering in RRTMG is ordered from surface to toa. 
!    c) calls CLDPROP_SW to set cloud optical depth based on input
!       cloud properties
!    d) calls SETCOEF_SW to calculate various quantities needed for 
!       the radiative transfer algorithm
!    e) calls SPCVRT to call the two-stream model that in turn 
!       calls TAUMOL to calculate gaseous optical depths for each 
!       of the 16 spectral bands and to perform the radiative transfer;
!    f) passes the calculated fluxes and cooling rates back to GCM
!
! Two modes of operation are possible:
!     The mode is chosen by using either rrtmg_sw.nomcica.f90 (to not use
!     McICA) or rrtmg_sw.f90 (to use McICA) to interface with a GCM.
!
!    1) Standard, single forward model calculation (imca = 0); this is 
!       valid only for clear sky or fully overcast clouds
!    2) Monte Carlo Independent Column Approximation (McICA, Pincus et al., 
!       JC, 2003) method is applied to the forward model calculation (imca = 1)
!       This method is valid for clear sky and full or partial cloud conditions.
!
! Two methods of cloud property input are possible:
!     Cloud properties can be input in one of two ways (controlled by input 
!     flags inflag, iceflag and liqflag; see text file rrtmg_sw_instructions
!     and subroutine rrtmg_sw_cldprop.f90 for further details):
!
!    1) Input cloud fraction, cloud optical depth, single scattering albedo 
!       and asymmetry parameter directly (inflgsw = 0)
!    2) Input cloud fraction and cloud physical properties: ice fracion,
!       ice and liquid particle sizes (inflgsw = 1 or 2);  
!       cloud optical properties are calculated by cldprop or cldprmc based
!       on input settings of iceflgsw and liqflgsw
!
! Two methods of aerosol property input are possible:
!     Aerosol properties can be input in one of two ways (controlled by input 
!     flag iaer, see text file rrtmg_sw_instructions for further details):
!
!    1) Input aerosol optical depth, single scattering albedo and asymmetry
!       parameter directly by layer and spectral band (iaer=10)
!    2) Input aerosol optical depth and 0.55 micron directly by layer and use
!       one or more of six ECMWF aerosol types (iaer=6)
!
!
! ------- Modifications -------
!
! This version of RRTMG_SW has been modified from RRTM_SW to use a reduced
! set of g-point intervals and a two-stream model for application to GCMs. 
!
!-- Original version (derived from RRTM_SW)
!     2002: AER. Inc.
!-- Conversion to F90 formatting; addition of 2-stream radiative transfer
!     Feb 2003: J.-J. Morcrette, ECMWF
!-- Additional modifications for GCM application
!     Aug 2003: M. J. Iacono, AER Inc.
!-- Total number of g-points reduced from 224 to 112.  Original
!   set of 224 can be restored by exchanging code in module parrrsw.f90 
!   and in file rrtmg_sw_init.f90.
!     Apr 2004: M. J. Iacono, AER, Inc.
!-- Modifications to include output for direct and diffuse 
!   downward fluxes.  There are output as "true" fluxes without
!   any delta scaling applied.  Code can be commented to exclude
!   this calculation in source file rrtmg_sw_spcvrt.f90.
!     Jan 2005: E. J. Mlawer, M. J. Iacono, AER, Inc.
!-- Reformatted for consistency with rrtmg_lw.
!     Feb 2007: M. J. Iacono, AER, Inc.
!-- Modifications to formatting to use assumed-shape arrays. 
!     Aug 2007: M. J. Iacono, AER, Inc.
!-- Modified to output direct and diffuse fluxes either with or without
!   delta scaling based on setting of idelm flag. 
!     Dec 2008: M. J. Iacono, AER, Inc.

! --------- Modules ---------

      use parrrsw, only : nbndsw, ngptsw, naerec, nstr, nmol, mxmol, &
                          jpband, jpb1, jpb2
      use rrsw_aer, only : rsrtaua, rsrpiza, rsrasya
      use rrsw_con, only : heatfac, oneminus, pi
      use rrsw_wvn, only : wavenum1, wavenum2

! ------- Declarations

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
      integer(kind=4), intent(inout) :: iaer         ! Aerosol option flag
                                                      !    0: No aerosol
                                                      !    6: ECMWF method
                                                      !    10:Input aerosol optical 
                                                      !       properties
      integer(kind=4), intent(inout) :: idelm        ! delta-m scaling flag
                                                      ! [0 = direct and diffuse fluxes are unscaled]
                                                      ! [1 = direct and diffuse fluxes are scaled]
                                                      ! (total downward fluxes are always delta scaled)

      real(kind=8), intent(in) :: play(ncol,nlay)          ! Layer pressures (hPa, mb)
                                                      !    Dimensions: (ncol,nlay)
      real(kind=8), intent(in) :: plev(ncol,nlay+1)          ! Interface pressures (hPa, mb)
                                                      !    Dimensions: (ncol,nlay+1)
      real(kind=8), intent(in) :: tlay(ncol,nlay)          ! Layer temperatures (K)
                                                      !    Dimensions: (ncol,nlay)
      real(kind=8), intent(in) :: tlev(ncol,nlay+1)          ! Interface temperatures (K)
                                                      !    Dimensions: (ncol,nlay+1)
      real(kind=8), intent(in) :: tsfc(ncol)            ! Surface temperature (K)
                                                      !    Dimensions: (ncol)
      real(kind=8), intent(in) :: h2ovmr(ncol,nlay)        ! H2O volume mixing ratio
                                                      !    Dimensions: (ncol,nlay)
      real(kind=8), intent(in) :: o3vmr(ncol,nlay)         ! O3 volume mixing ratio
                                                      !    Dimensions: (ncol,nlay)
      real(kind=8), intent(in) :: co2vmr(ncol,nlay)        ! CO2 volume mixing ratio
                                                      !    Dimensions: (ncol,nlay)
      real(kind=8), intent(in) :: ch4vmr(ncol,nlay)        ! Methane volume mixing ratio
                                                      !    Dimensions: (ncol,nlay)
      real(kind=8), intent(in) :: n2ovmr(ncol,nlay)        ! Nitrous oxide volume mixing ratio
                                                      !    Dimensions: (ncol,nlay)
      real(kind=8), intent(in) :: o2vmr(ncol,nlay)         ! Oxygen volume mixing ratio
                                                      !    Dimensions: (ncol,nlay)
      real(kind=8), intent(in) :: asdir(ncol)           ! UV/vis surface albedo direct rad
                                                      !    Dimensions: (ncol)
      real(kind=8), intent(in) :: aldir(ncol)           ! Near-IR surface albedo direct rad
                                                      !    Dimensions: (ncol)
      real(kind=8), intent(in) :: asdif(ncol)           ! UV/vis surface albedo: diffuse rad
                                                      !    Dimensions: (ncol)
      real(kind=8), intent(in) :: aldif(ncol)           ! Near-IR surface albedo: diffuse rad
                                                      !    Dimensions: (ncol)

      integer(kind=4), intent(in) :: dyofyr          ! Day of the year (used to get Earth/Sun
                                                      !  distance if adjflx not provided)
      real(kind=8), intent(in) :: adjes              ! Flux adjustment for Earth/Sun distance
      real(kind=8), intent(in) :: coszen(ncol)          ! Cosine of solar zenith angle
                                                      !    Dimensions: (ncol)
      real(kind=8), intent(in) :: scon               ! Solar constant (W/m2)

      integer(kind=4), intent(in) :: inflgsw         ! Flag for cloud optical properties
      integer(kind=4), intent(in) :: iceflgsw        ! Flag for ice particle specification
      integer(kind=4), intent(in) :: liqflgsw        ! Flag for liquid droplet specification

      real(kind=8), intent(in) :: cldfr(ncol,nlay)         ! Cloud fraction
                                                      !    Dimensions: (ncol,nlay)
      real(kind=8), intent(in) :: taucld(nbndsw,ncol,nlay)      ! In-cloud optical depth
                                                      !    Dimensions: (nbndsw,ncol,nlay)
      real(kind=8), intent(in) :: ssacld(nbndsw,ncol,nlay)      ! In-cloud single scattering albedo
                                                      !    Dimensions: (nbndsw,ncol,nlay)
      real(kind=8), intent(in) :: asmcld(nbndsw,ncol,nlay)      ! In-cloud asymmetry parameter
                                                      !    Dimensions: (nbndsw,ncol,nlay)
      real(kind=8), intent(in) :: fsfcld(nbndsw,ncol,nlay)      ! In-cloud forward scattering fraction
                                                      !    Dimensions: (nbndsw,ncol,nlay)
      real(kind=8), intent(in) :: cicewp(ncol,nlay)        ! In-cloud ice water path (g/m2)
                                                      !    Dimensions: (ncol,nlay)
      real(kind=8), intent(in) :: cliqwp(ncol,nlay)        ! In-cloud liquid water path (g/m2)
                                                      !    Dimensions: (ncol,nlay)
      real(kind=8), intent(in) :: reice(ncol,nlay)         ! Cloud ice effective radius (microns)
                                                      !    Dimensions: (ncol,nlay)
                                                      ! specific definition of reice depends on setting of iceflgsw:
                                                      ! iceflgsw = 0: (inactive)
                                                      !              
                                                      ! iceflgsw = 1: ice effective radius, r_ec, (Ebert and Curry, 1992),
                                                      !               r_ec range is limited to 13.0 to 130.0 microns
                                                      ! iceflgsw = 2: ice effective radius, r_k, (Key, Streamer Ref. Manual, 1996)
                                                      !               r_k range is limited to 5.0 to 131.0 microns
                                                      ! iceflgsw = 3: generalized effective size, dge, (Fu, 1996),
                                                      !               dge range is limited to 5.0 to 140.0 microns
                                                      !               [dge = 1.0315 * r_ec]
      real(kind=8), intent(in) :: reliq(ncol,nlay)         ! Cloud water drop effective radius (microns)
                                                      !    Dimensions: (ncol,nlay)
      real(kind=8), intent(in) :: tauaer(ncol,nlay,nbndsw)      ! Aerosol optical depth (iaer=10 only)
                                                      !    Dimensions: (ncol,nlay,nbndsw)
                                                      ! (non-delta scaled)      
      real(kind=8), intent(in) :: ssaaer(ncol,nlay,nbndsw)      ! Aerosol single scattering albedo (iaer=10 only)
                                                      !    Dimensions: (ncol,nlay,nbndsw)
                                                      ! (non-delta scaled)      
      real(kind=8), intent(in) :: asmaer(ncol,nlay,nbndsw)      ! Aerosol asymmetry parameter (iaer=10 only)
                                                      !    Dimensions: (ncol,nlay,nbndsw)
                                                      ! (non-delta scaled)      
      real(kind=8), intent(in) :: ecaer(ncol,nlay,naerec)       ! Aerosol optical depth at 0.55 micron (iaer=6 only)
                                                      !    Dimensions: (ncol,nlay,naerec)
                                                      ! (non-delta scaled)      

! ----- Output -----

      real(kind=8), intent(out) :: swuflx(ncol,nlay+1)       ! Total sky shortwave upward flux (W/m2)
                                                      !    Dimensions: (ncol,nlay+1)
      real(kind=8), intent(out) :: swdflx(ncol,nlay+1)       ! Total sky shortwave downward flux (W/m2)
                                                      !    Dimensions: (ncol,nlay+1)
      real(kind=8), intent(out) :: swhr(ncol,nlay)         ! Total sky shortwave radiative heating rate (K/d)
                                                      !    Dimensions: (ncol,nlay)
      real(kind=8), intent(out) :: swuflxc(ncol,nlay+1)      ! Clear sky shortwave upward flux (W/m2)
                                                      !    Dimensions: (ncol,nlay+1)
      real(kind=8), intent(out) :: swdflxc(ncol,nlay+1)      ! Clear sky shortwave downward flux (W/m2)
                                                      !    Dimensions: (ncol,nlay+1)
      real(kind=8), intent(out) :: swhrc(ncol,nlay)        ! Clear sky shortwave radiative heating rate (K/d)
                                                      !    Dimensions: (ncol,nlay)
      ! Added Output by ZTAN
      real(kind=8), intent(out) :: dirdflux(ncol,nlay+1)       ! Direct downward shortwave flux (W/m2)
                                                      !    Dimensions: (ncol,nlay+1)
      real(kind=8), intent(out) :: difdflux(ncol,nlay+1)       ! Diffuse downward shortwave flux (W/m2)
                                                      !    Dimensions: (ncol,nlay+1)
                                                      
! ----- Local -----

! Control
      integer(kind=4) :: nlayers             ! total number of layers
      integer(kind=4) :: istart              ! beginning band of calculation
      integer(kind=4) :: iend                ! ending band of calculation
      integer(kind=4) :: icpr                ! cldprop/cldprmc use flag
      integer(kind=4) :: iout                ! output option flag
!      integer(kind=4) :: idelm               ! delta-m scaling flag
!                                              ! [0 = direct and diffuse fluxes are unscaled]
!                                              ! [1 = direct and diffuse fluxes are scaled]
!                                              ! (total downward fluxes are always delta scaled)
      integer(kind=4) :: isccos              ! instrumental cosine response flag (inactive)
      integer(kind=4) :: iplon               ! column loop index
      integer(kind=4) :: i                   ! layer loop index                       ! jk
      integer(kind=4) :: ib                  ! band loop index                        ! jsw
      integer(kind=4) :: ia, ig              ! indices
      integer(kind=4) :: k                   ! layer loop index
      integer(kind=4) :: ims                 ! value for changing mcica permute seed
      integer(kind=4) :: imca                ! flag for mcica [0=off, 1=on]

      real(kind=8) :: zepsec, zepzen         ! epsilon
      real(kind=8) :: zdpgcp                 ! flux to heating conversion ratio

! Atmosphere
      real(kind=8) :: pavel(nlay)          ! layer pressures (mb) 
      real(kind=8) :: tavel(nlay)          ! layer temperatures (K)
      real(kind=8) :: pz(0:nlay)           ! level (interface) pressures (hPa, mb)
      real(kind=8) :: tz(0:nlay)           ! level (interface) temperatures (K)
      real(kind=8) :: tbound                 ! surface temperature (K)
      real(kind=8) :: pdp(nlay)            ! layer pressure thickness (hPa, mb)
      real(kind=8) :: coldry(nlay)         ! dry air column amount
      real(kind=8) :: wkl(mxmol,nlay)      ! molecular amounts (mol/cm-2)

!      real(kind=8) :: earth_sun             ! function for Earth/Sun distance factor
      real(kind=8) :: cossza                 ! Cosine of solar zenith angle
      real(kind=8) :: adjflux(jpband)        ! adjustment for current Earth/Sun distance
      real(kind=8) :: solvar(jpband)         ! solar constant scaling factor from rrtmg_sw
                                              !  default value of 1368.22 Wm-2 at 1 AU
      real(kind=8) :: albdir(nbndsw)         ! surface albedo, direct          ! zalbp
      real(kind=8) :: albdif(nbndsw)         ! surface albedo, diffuse         ! zalbd

      real(kind=8) :: taua(nlay,nbndsw)    ! Aerosol optical depth
      real(kind=8) :: ssaa(nlay,nbndsw)    ! Aerosol single scattering albedo
      real(kind=8) :: asma(nlay,nbndsw)    ! Aerosol asymmetry parameter

! Atmosphere - setcoef
      integer(kind=4) :: laytrop             ! tropopause layer index
      integer(kind=4) :: layswtch            ! tropopause layer index
      integer(kind=4) :: laylow              ! tropopause layer index
      integer(kind=4) :: jp(nlay)          ! 
      integer(kind=4) :: jt(nlay)          !
      integer(kind=4) :: jt1(nlay)         !

      real(kind=8) :: colh2o(nlay)         ! column amount (h2o)
      real(kind=8) :: colco2(nlay)         ! column amount (co2)
      real(kind=8) :: colo3(nlay)          ! column amount (o3)
      real(kind=8) :: coln2o(nlay)         ! column amount (n2o)
      real(kind=8) :: colch4(nlay)         ! column amount (ch4)
      real(kind=8) :: colo2(nlay)          ! column amount (o2)
      real(kind=8) :: colmol(nlay)         ! column amount
      real(kind=8) :: co2mult(nlay)        ! column amount 

      integer(kind=4) :: indself(nlay)
      integer(kind=4) :: indfor(nlay)
      real(kind=8) :: selffac(nlay)
      real(kind=8) :: selffrac(nlay)
      real(kind=8) :: forfac(nlay)
      real(kind=8) :: forfrac(nlay)

      real(kind=8) :: &                      !
                         fac00(nlay), fac01(nlay), &
                         fac10(nlay), fac11(nlay) 

! Atmosphere/clouds - cldprop
      integer(kind=4) :: ncbands             ! number of cloud spectral bands
      integer(kind=4) :: inflag              ! flag for cloud property method
      integer(kind=4) :: iceflag             ! flag for ice cloud properties
      integer(kind=4) :: liqflag             ! flag for liquid cloud properties

      real(kind=8) :: cldfrac(nlay)        ! layer cloud fraction
      real(kind=8) :: tauc(nbndsw,nlay)    ! in-cloud optical depth (non-delta scaled)
      real(kind=8) :: ssac(nbndsw,nlay)    ! in-cloud single scattering albedo (non-delta scaled)
      real(kind=8) :: asmc(nbndsw,nlay)    ! in-cloud asymmetry parameter (non-delta scaled)
      real(kind=8) :: fsfc(nbndsw,nlay)    ! in-cloud forward scattering fraction (non-delta scaled)
      real(kind=8) :: ciwp(nlay)           ! in-cloud ice water path
      real(kind=8) :: clwp(nlay)           ! in-cloud liquid water path
      real(kind=8) :: rel(nlay)            ! cloud liquid particle effective radius (microns)
      real(kind=8) :: rei(nlay)            ! cloud ice particle effective size (microns)

      real(kind=8) :: taucloud(nlay,jpband)  ! in-cloud optical depth
      real(kind=8) :: taucldorig(nlay,jpband)! in-cloud optical depth (non-delta scaled)
      real(kind=8) :: ssacloud(nlay,jpband)  ! in-cloud single scattering albedo
      real(kind=8) :: asmcloud(nlay,jpband)  ! in-cloud asymmetry parameter

! Atmosphere/clouds/aerosol - spcvrt,spcvmc
      real(kind=8) :: ztauc(nlay,nbndsw)     ! cloud optical depth
      real(kind=8) :: ztaucorig(nlay,nbndsw) ! unscaled cloud optical depth
      real(kind=8) :: zasyc(nlay,nbndsw)     ! cloud asymmetry parameter 
                                                !  (first moment of phase function)
      real(kind=8) :: zomgc(nlay,nbndsw)     ! cloud single scattering albedo
      real(kind=8) :: ztaua(nlay,nbndsw)     ! total aerosol optical depth
      real(kind=8) :: zasya(nlay,nbndsw)     ! total aerosol asymmetry parameter 
      real(kind=8) :: zomga(nlay,nbndsw)     ! total aerosol single scattering albedo

      real(kind=8) :: zbbfu(nlay+1)          ! temporary upward shortwave flux (w/m2)
      real(kind=8) :: zbbfd(nlay+1)          ! temporary downward shortwave flux (w/m2)
      real(kind=8) :: zbbcu(nlay+1)          ! temporary clear sky upward shortwave flux (w/m2)
      real(kind=8) :: zbbcd(nlay+1)          ! temporary clear sky downward shortwave flux (w/m2)
      real(kind=8) :: zbbfddir(nlay+1)       ! temporary downward direct shortwave flux (w/m2)
      real(kind=8) :: zbbcddir(nlay+1)       ! temporary clear sky downward direct shortwave flux (w/m2)
      real(kind=8) :: zuvfd(nlay+1)          ! temporary UV downward shortwave flux (w/m2)
      real(kind=8) :: zuvcd(nlay+1)          ! temporary clear sky UV downward shortwave flux (w/m2)
      real(kind=8) :: zuvfddir(nlay+1)       ! temporary UV downward direct shortwave flux (w/m2)
      real(kind=8) :: zuvcddir(nlay+1)       ! temporary clear sky UV downward direct shortwave flux (w/m2)
      real(kind=8) :: znifd(nlay+1)          ! temporary near-IR downward shortwave flux (w/m2)
      real(kind=8) :: znicd(nlay+1)          ! temporary clear sky near-IR downward shortwave flux (w/m2)
      real(kind=8) :: znifddir(nlay+1)       ! temporary near-IR downward direct shortwave flux (w/m2)
      real(kind=8) :: znicddir(nlay+1)       ! temporary clear sky near-IR downward direct shortwave flux (w/m2)

! Optional output fields 
      real(kind=8) :: swnflx(nlay+1)         ! Total sky shortwave net flux (W/m2)
      real(kind=8) :: swnflxc(nlay+1)        ! Clear sky shortwave net flux (W/m2)
      ! real(kind=8) :: dirdflux(nlay+1)       ! Direct downward shortwave surface flux
      ! real(kind=8) :: difdflux(nlay+1)       ! Diffuse downward shortwave surface flux
      real(kind=8) :: uvdflx(nlay+1)         ! Total sky downward shortwave flux, UV/vis  
      real(kind=8) :: nidflx(nlay+1)         ! Total sky downward shortwave flux, near-IR 
      real(kind=8) :: dirdnuv(nlay+1)        ! Direct downward shortwave flux, UV/vis
      real(kind=8) :: difdnuv(nlay+1)        ! Diffuse downward shortwave flux, UV/vis
      real(kind=8) :: dirdnir(nlay+1)        ! Direct downward shortwave flux, near-IR
      real(kind=8) :: difdnir(nlay+1)        ! Diffuse downward shortwave flux, near-IR

! Output - inactive
!      real(kind=8) :: zuvfu(nlay+1)         ! temporary upward UV shortwave flux (w/m2)
!      real(kind=8) :: zuvfd(nlay+1)         ! temporary downward UV shortwave flux (w/m2)
!      real(kind=8) :: zuvcu(nlay+1)         ! temporary clear sky upward UV shortwave flux (w/m2)
!      real(kind=8) :: zuvcd(nlay+1)         ! temporary clear sky downward UV shortwave flux (w/m2)
!      real(kind=8) :: zvsfu(nlay+1)         ! temporary upward visible shortwave flux (w/m2)
!      real(kind=8) :: zvsfd(nlay+1)         ! temporary downward visible shortwave flux (w/m2)
!      real(kind=8) :: zvscu(nlay+1)         ! temporary clear sky upward visible shortwave flux (w/m2)
!      real(kind=8) :: zvscd(nlay+1)         ! temporary clear sky downward visible shortwave flux (w/m2)
!      real(kind=8) :: znifu(nlay+1)         ! temporary upward near-IR shortwave flux (w/m2)
!      real(kind=8) :: znifd(nlay+1)         ! temporary downward near-IR shortwave flux (w/m2)
!      real(kind=8) :: znicu(nlay+1)         ! temporary clear sky upward near-IR shortwave flux (w/m2)
!      real(kind=8) :: znicd(nlay+1)         ! temporary clear sky downward near-IR shortwave flux (w/m2)


! Initializations

      zepsec = 1.e-06_8
      zepzen = 1.e-10_8
      oneminus = 1.0_8 - zepsec
      pi = 2._8 * asin(1._8)

      istart = jpb1
      iend = jpb2
      iout = 0
      icpr = 0

! In a GCM with or without McICA, set nlon to the longitude dimension
!
! Set imca to select calculation type:
!  imca = 0, use standard forward model calculation (clear and overcast only)
!  imca = 1, use McICA for Monte Carlo treatment of sub-grid cloud variability
!            (clear, overcast or partial cloud conditions)

! *** This version does not use McICA (imca = 0) ***

! Set icld to select of clear or cloud calculation and cloud 
! overlap method (read by subroutine readprof from input file INPUT_RRTM):  
! Without McICA, SW calculation is limited to clear or fully overcast conditions. 
! icld = 0, clear only
! icld = 1, with clouds using random cloud overlap (McICA only)
! icld = 2, with clouds using maximum/random cloud overlap (McICA only)
! icld = 3, with clouds using maximum cloud overlap (McICA only)
      if (icld.lt.0.or.icld.gt.3) icld = 2

! Set iaer to select aerosol option
! iaer = 0, no aerosols
! iaer = 6, use six ECMWF aerosol types
!           input aerosol optical depth at 0.55 microns for each aerosol type (ecaer)
! iaer = 10, input total aerosol optical depth, single scattering albedo 
!            and asymmetry parameter (tauaer, ssaaer, asmaer) directly
      if (iaer.ne.0.and.iaer.ne.6.and.iaer.ne.10) iaer = 0

! Set idelm to select between delta-M scaled or unscaled output direct and diffuse fluxes
! NOTE: total downward fluxes are always delta scaled
! idelm = 0, output direct and diffuse flux components are not delta scaled
!            (direct flux does not include forward scattering peak)
! idelm = 1, output direct and diffuse flux components are delta scaled (default)
!            (direct flux includes part or most of forward scattering peak)
!      idelm = 1
     if (idelm.ne.0.and.idelm.ne.1)  idelm = 1

! Call model and data initialization, compute lookup tables, perform
! reduction of g-points from 224 to 112 for input absorption
! coefficient data and other arrays.
!
! In a GCM this call should be placed in the model initialization
! area, since this has to be called only once.  
!      call rrtmg_sw_ini(cpdair)

! This is the main longitude/column loop in RRTMG.
! Modify to loop over all columns (nlon) or over daylight columns

      do iplon = 1, ncol

! Prepare atmosphere profile from GCM for use in RRTMG, and define
! other input parameters

         call inatm_sw (iplon, ncol, nlay, icld, iaer, &
              play, plev, tlay, tlev, tsfc, h2ovmr, &
              o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr, &
              adjes, dyofyr, scon, inflgsw, iceflgsw, liqflgsw, &
              cldfr, taucld, ssacld, asmcld, fsfcld, cicewp, cliqwp, &
              reice, reliq, tauaer, ssaaer, asmaer, &
              nlayers, pavel, pz, pdp, tavel, tz, tbound, coldry, wkl, &
              adjflux, solvar, inflag, iceflag, liqflag, cldfrac, tauc, &
              ssac, asmc, fsfc, ciwp, clwp, rei, rel, taua, ssaa, asma)

!  For cloudy atmosphere, use cldprop to set cloud optical properties based on
!  input cloud physical properties.  Select method based on choices described
!  in cldprop.  Cloud fraction, water path, liquid droplet and ice particle
!  effective radius must be passed in cldprop.  Cloud fraction and cloud
!  optical properties are transferred to rrtmg_sw arrays in cldprop.  

!  Without McICA, SW calculation is limited to clear or fully overcast conditions. 
!  Stop model if partial cloudiness is present.

         do i = 1, nlayers
            if (cldfrac(i).gt.zepsec .and. cldfrac(i).lt.oneminus) then
               stop 'PARTIAL CLOUD NOT ALLOWED'
            endif
         enddo
         call cldprop_sw(nlayers, inflag, iceflag, liqflag, cldfrac, &
                         tauc, ssac, asmc, fsfc, ciwp, clwp, rei, rel, &
                         taucldorig, taucloud, ssacloud, asmcloud)
         icpr = 1  

! Calculate coefficients for the temperature and pressure dependence of the 
! molecular absorption coefficients by interpolating data from stored
! reference atmospheres.

         call setcoef_sw(nlayers, pavel, tavel, pz, tz, tbound, coldry, wkl, &
                         laytrop, layswtch, laylow, jp, jt, jt1, &
                         co2mult, colch4, colco2, colh2o, colmol, coln2o, &
                         colo2, colo3, fac00, fac01, fac10, fac11, &
                         selffac, selffrac, indself, forfac, forfrac, indfor)


! Cosine of the solar zenith angle 
!  Prevent using value of zero; ideally, SW model is not called from host model when sun 
!  is below horizon

         cossza = coszen(iplon)
         if (cossza .lt. zepzen) cossza = zepzen


! Transfer albedo, cloud and aerosol properties into arrays for 2-stream radiative transfer 

! Surface albedo
!  Near-IR bands 16-24 and 29 (1-9 and 14), 820-16000 cm-1, 0.625-12.195 microns
         do ib=1,9
            albdir(ib) = aldir(iplon)
            albdif(ib) = aldif(iplon)
         enddo
         albdir(nbndsw) = aldir(iplon)
         albdif(nbndsw) = aldif(iplon)
!  UV/visible bands 25-28 (10-13), 16000-50000 cm-1, 0.200-0.625 micron
         do ib=10,13
            albdir(ib) = asdir(iplon)
            albdif(ib) = asdif(iplon)
         enddo


! Clouds
         if (icld.eq.0) then

            ztauc(:,:) = 0._8
            ztaucorig(:,:) = 0._8
            zasyc(:,:) = 0._8
            zomgc(:,:) = 1._8

         elseif (icld.ge.1) then
            do i=1,nlayers
               do ib=1,nbndsw
                  ztauc(i,ib) = taucloud(i,jpb1-1+ib)
                  ztaucorig(i,ib) = taucldorig(i,jpb1-1+ib)
                  zasyc(i,ib) = asmcloud(i,jpb1-1+ib)
                  zomgc(i,ib) = ssacloud(i,jpb1-1+ib)
               enddo
            enddo

         endif   

! Aerosol
! IAER = 0: no aerosols
         if (iaer.eq.0) then

            ztaua(:,:) = 0._8
            zasya(:,:) = 0._8
            zomga(:,:) = 1._8

! IAER = 6: Use ECMWF six aerosol types. See rrsw_aer.f90 for details.
! Input aerosol optical thickness at 0.55 micron for each aerosol type (ecaer), 
! or set manually here for each aerosol and layer.
         elseif (iaer.eq.6) then

!            do i = 1, nlayers
!               do ia = 1, naerec
!                  ecaer(iplon,i,ia) = 1.0e-15_8
!               enddo
!            enddo

            do i = 1, nlayers
               do ib = 1, nbndsw
                  ztaua(i,ib) = 0._8
                  zasya(i,ib) = 0._8
                  zomga(i,ib) = 0._8
                  do ia = 1, naerec
                     ztaua(i,ib) = ztaua(i,ib) + rsrtaua(ib,ia) * ecaer(iplon,i,ia)
                     zomga(i,ib) = zomga(i,ib) + rsrtaua(ib,ia) * ecaer(iplon,i,ia) * &
                                   rsrpiza(ib,ia)
                     zasya(i,ib) = zasya(i,ib) + rsrtaua(ib,ia) * ecaer(iplon,i,ia) * &
                                   rsrpiza(ib,ia) * rsrasya(ib,ia)
                  enddo
                  if (ztaua(i,ib) == 0._8) then
                     ztaua(i,ib) = 0._8
                     zasya(i,ib) = 0._8
                     zomga(i,ib) = 1._8
                  else
                     if (zomga(i,ib) /= 0._8) then
                        zasya(i,ib) = zasya(i,ib) / zomga(i,ib)
                     endif
                     if (ztaua(i,ib) /= 0._8) then
                        zomga(i,ib) = zomga(i,ib) / ztaua(i,ib)
                     endif
                  endif
               enddo
            enddo

! IAER=10: Direct specification of aerosol optical properties from GCM
         elseif (iaer.eq.10) then

            do i = 1 ,nlayers
               do ib = 1 ,nbndsw
                  ztaua(i,ib) = taua(i,ib)
                  zasya(i,ib) = asma(i,ib)
                  zomga(i,ib) = ssaa(i,ib)
               enddo
            enddo

         endif


! Call the 2-stream radiation transfer model

         do i=1,nlayers+1
            zbbcu(i) = 0._8
            zbbcd(i) = 0._8
            zbbfu(i) = 0._8
            zbbfd(i) = 0._8
            zbbcddir(i) = 0._8
            zbbfddir(i) = 0._8
            zuvcd(i) = 0._8
            zuvfd(i) = 0._8
            zuvcddir(i) = 0._8
            zuvfddir(i) = 0._8
            znicd(i) = 0._8
            znifd(i) = 0._8
            znicddir(i) = 0._8
            znifddir(i) = 0._8
         enddo

         ! write(*,*) ztaua(1:12, 1:nbndsw)
         ! write(*,*) zasya(1:12, 1:nbndsw)
         ! write(*,*) zomga(1:12, 1:nbndsw)
            
         call spcvrt_sw &
             (nlayers, istart, iend, icpr, idelm, iout, &
              pavel, tavel, pz, tz, tbound, albdif, albdir, &
              cldfrac, ztauc, zasyc, zomgc, ztaucorig, &
              ztaua, zasya, zomga, cossza, coldry, wkl, adjflux, &	 
              laytrop, layswtch, laylow, jp, jt, jt1, &
              co2mult, colch4, colco2, colh2o, colmol, coln2o, colo2, colo3, &
              fac00, fac01, fac10, fac11, &
              selffac, selffrac, indself, forfac, forfrac, indfor, &
              zbbfd, zbbfu, zbbcd, zbbcu, zuvfd, zuvcd, znifd, znicd, &
              zbbfddir, zbbcddir, zuvfddir, zuvcddir, znifddir, znicddir)
          ! write(*,*) zbbfddir

! Transfer up and down, clear and total sky fluxes to output arrays.
! Vertical indexing goes from bottom to top; reverse here for GCM if necessary.

         do i = 1, nlayers+1
            swuflxc(iplon,i) = zbbcu(i)
            swdflxc(iplon,i) = zbbcd(i)
            swuflx(iplon,i) = zbbfu(i)
            swdflx(iplon,i) = zbbfd(i)
            uvdflx(i) = zuvfd(i)
            nidflx(i) = znifd(i)
!  Direct/diffuse fluxes
            ! dirdflux(i) = zbbfddir(i)
            ! difdflux(i) = swdflx(iplon,i) - dirdflux(i)
            dirdflux(iplon,i) = zbbfddir(i)
            difdflux(iplon,i) = swdflx(iplon,i) - dirdflux(iplon,i)
!  UV/visible direct/diffuse fluxes
            dirdnuv(i) = zuvfddir(i)
            difdnuv(i) = zuvfd(i) - dirdnuv(i)
!  Near-IR direct/diffuse fluxes
            dirdnir(i) = znifddir(i)
            difdnir(i) = znifd(i) - dirdnir(i)
         enddo

!  Total and clear sky net fluxes
         do i = 1, nlayers+1
            swnflxc(i) = swdflxc(iplon,i) - swuflxc(iplon,i)
            swnflx(i) = swdflx(iplon,i) - swuflx(iplon,i)
         enddo

!  Total and clear sky heating rates
         do i = 1, nlayers
            zdpgcp = heatfac / pdp(i)
            swhrc(iplon,i) = (swnflxc(i+1) - swnflxc(i)) * zdpgcp
            swhr(iplon,i) = (swnflx(i+1) - swnflx(i)) * zdpgcp
         enddo
         !swhrc(iplon,nlayers) = 0._8
         !swhr(iplon,nlayers) = 0._8

! End longitude loop
      enddo

      end subroutine rrtmg_sw

!*************************************************************************
      real(kind=8) function earth_sun(idn)
!*************************************************************************
!
!  Purpose: Function to calculate the correction factor of Earth's orbit
!  for current day of the year

!  idn        : Day of the year
!  earth_sun  : square of the ratio of mean to actual Earth-Sun distance

! ------- Modules -------

      use rrsw_con, only : pi

      integer(kind=4), intent(in) :: idn

      real(kind=8) :: gamma

      gamma = 2._8*pi*(idn-1)/365._8

! Use Iqbal's equation 1.2.1

      earth_sun = 1.000110_8 + .034221_8 * cos(gamma) + .001289_8 * sin(gamma) + &
                   .000719_8 * cos(2._8*gamma) + .000077_8 * sin(2._8*gamma)

      end function earth_sun

!***************************************************************************
      subroutine inatm_sw (iplon, ncol, nlay, icld, iaer, &
            play, plev, tlay, tlev, tsfc, h2ovmr, &
            o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr, &
            adjes, dyofyr, scon, inflgsw, iceflgsw, liqflgsw, &
            cldfr, taucld, ssacld, asmcld, fsfcld, cicewp, cliqwp, &
            reice, reliq, tauaer, ssaaer, asmaer, &
            nlayers, pavel, pz, pdp, tavel, tz, tbound, coldry, wkl, &
            adjflux, solvar, inflag, iceflag, liqflag, cldfrac, tauc, &
            ssac, asmc, fsfc, ciwp, clwp, rei, rel, taua, ssaa, asma)
!***************************************************************************
!
!  Input atmospheric profile from GCM, and prepare it for use in RRTMG_SW.
!  Set other RRTMG_SW input parameters.  
!
!***************************************************************************

! --------- Modules ----------

      use parrrsw, only : nbndsw, ngptsw, nstr, nmol, mxmol, &
                          jpband, jpb1, jpb2, rrsw_scon
      use rrsw_con, only : fluxfac, heatfac, oneminus, pi, grav, avogad
      use rrsw_wvn, only : ng, nspa, nspb, wavenum1, wavenum2, delwave

! ------- Declarations -------

! ----- Input -----
! Note: All volume mixing ratios are in dimensionless units of mole fraction obtained
! by scaling mass mixing ratio (g/g) with the appropriate molecular weights (g/mol) 
      integer(kind=4), intent(in) :: iplon           ! column loop index
      integer(kind=4), intent(in) :: ncol            ! Number of horizontal columns
      integer(kind=4), intent(in) :: nlay            ! number of model layers
      integer(kind=4), intent(in) :: icld            ! clear/cloud flag
      integer(kind=4), intent(in) :: iaer            ! aerosol option flag

      real(kind=8), intent(in) :: play(ncol,nlay)          ! Layer pressures (hPa, mb)
                                                      ! Dimensions: (ncol,nlay)
      real(kind=8), intent(in) :: plev(ncol,nlay+1)          ! Interface pressures (hPa, mb)
                                                      ! Dimensions: (ncol,nlay+1)
      real(kind=8), intent(in) :: tlay(ncol,nlay)          ! Layer temperatures (K)
                                                      ! Dimensions: (ncol,nlay)
      real(kind=8), intent(in) :: tlev(ncol,nlay+1)          ! Interface temperatures (K)
                                                      ! Dimensions: (ncol,nlay+1)
      real(kind=8), intent(in) :: tsfc(ncol)            ! Surface temperature (K)
                                                      ! Dimensions: (ncol)
      real(kind=8), intent(in) :: h2ovmr(ncol,nlay)        ! H2O volume mixing ratio
                                                      ! Dimensions: (ncol,nlay)
      real(kind=8), intent(in) :: o3vmr(ncol,nlay)         ! O3 volume mixing ratio
                                                      ! Dimensions: (ncol,nlay)
      real(kind=8), intent(in) :: co2vmr(ncol,nlay)        ! CO2 volume mixing ratio
                                                      ! Dimensions: (ncol,nlay)
      real(kind=8), intent(in) :: ch4vmr(ncol,nlay)        ! Methane volume mixing ratio
                                                      ! Dimensions: (ncol,nlay)
      real(kind=8), intent(in) :: n2ovmr(ncol,nlay)        ! Nitrous oxide volume mixing ratio
                                                      ! Dimensions: (ncol,nlay)
      real(kind=8), intent(in) :: o2vmr(ncol,nlay)         ! Oxygen volume mixing ratio
                                                      ! Dimensions: (ncol,nlay)

      integer(kind=4), intent(in) :: dyofyr          ! Day of the year (used to get Earth/Sun
                                                      !  distance if adjflx not provided)
      real(kind=8), intent(in) :: adjes              ! Flux adjustment for Earth/Sun distance
      real(kind=8), intent(in) :: scon               ! Solar constant (W/m2)

      integer(kind=4), intent(in) :: inflgsw         ! Flag for cloud optical properties
      integer(kind=4), intent(in) :: iceflgsw        ! Flag for ice particle specification
      integer(kind=4), intent(in) :: liqflgsw        ! Flag for liquid droplet specification

      real(kind=8), intent(in) :: cldfr(ncol,nlay)         ! Cloud fraction
                                                      ! Dimensions: (ncol,nlay)
      real(kind=8), intent(in) :: taucld(nbndsw,ncol,nlay)      ! In-cloud optical depth (optional)
                                                      ! Dimensions: (nbndsw,ncol,nlay)
      real(kind=8), intent(in) :: ssacld(nbndsw,ncol,nlay)      ! In-cloud single scattering albedo
                                                      ! Dimensions: (nbndsw,ncol,nlay)
      real(kind=8), intent(in) :: asmcld(nbndsw,ncol,nlay)      ! In-cloud asymmetry parameter
                                                      ! Dimensions: (nbndsw,ncol,nlay)
      real(kind=8), intent(in) :: fsfcld(nbndsw,ncol,nlay)      ! In-cloud forward scattering fraction
                                                      ! Dimensions: (nbndsw,ncol,nlay)
      real(kind=8), intent(in) :: cicewp(ncol,nlay)        ! In-cloud ice water path (g/m2)
                                                      ! Dimensions: (ncol,nlay)
      real(kind=8), intent(in) :: cliqwp(ncol,nlay)        ! In-cloud liquid water path (g/m2)
                                                      ! Dimensions: (ncol,nlay)
      real(kind=8), intent(in) :: reice(ncol,nlay)         ! Cloud ice effective size (microns)
                                                      ! Dimensions: (ncol,nlay)
      real(kind=8), intent(in) :: reliq(ncol,nlay)         ! Cloud water drop effective radius (microns)
                                                      ! Dimensions: (ncol,nlay)

      real(kind=8), intent(in) :: tauaer(ncol,nlay,nbndsw)      ! Aerosol optical depth
                                                      ! Dimensions: (ncol,nlay,nbndsw)
      real(kind=8), intent(in) :: ssaaer(ncol,nlay,nbndsw)      ! Aerosol single scattering albedo
                                                      ! Dimensions: (ncol,nlay,nbndsw)
      real(kind=8), intent(in) :: asmaer(ncol,nlay,nbndsw)      ! Aerosol asymmetry parameter
                                                      ! Dimensions: (ncol,nlay,nbndsw)

! Atmosphere
      integer(kind=4), intent(out) :: nlayers        ! number of layers

      real(kind=8), intent(out) :: pavel(nlay)          ! layer pressures (mb) 
                                                      ! Dimensions: (nlay)
      real(kind=8), intent(out) :: tavel(nlay)          ! layer temperatures (K)
                                                      ! Dimensions: (nlay)
      real(kind=8), intent(out) :: pz(0:nlay)            ! level (interface) pressures (hPa, mb)
                                                      ! Dimensions: (0:nlay)
      real(kind=8), intent(out) :: tz(0:nlay)            ! level (interface) temperatures (K)
                                                      ! Dimensions: (0:nlay)
      real(kind=8), intent(out) :: tbound            ! surface temperature (K)
      real(kind=8), intent(out) :: pdp(nlay)            ! layer pressure thickness (hPa, mb)
                                                      ! Dimensions: (nlay)
      real(kind=8), intent(out) :: coldry(nlay)         ! dry air column density (mol/cm2)
                                                      ! Dimensions: (nlay)
      real(kind=8), intent(out) :: wkl(mxmol,nlay)          ! molecular amounts (mol/cm-2)
                                                      ! Dimensions: (mxmol,nlay)

      real(kind=8), intent(out) :: adjflux(jpband)        ! adjustment for current Earth/Sun distance
                                                      ! Dimensions: (jpband)
      real(kind=8), intent(out) :: solvar(jpband)         ! solar constant scaling factor from rrtmg_sw
                                                      ! Dimensions: (jpband)
                                                      !  default value of 1368.22 Wm-2 at 1 AU
      real(kind=8), intent(out) :: taua(nlay,nbndsw)         ! Aerosol optical depth
                                                      ! Dimensions: (nlay,nbndsw)
      real(kind=8), intent(out) :: ssaa(nlay,nbndsw)         ! Aerosol single scattering albedo
                                                      ! Dimensions: (nlay,nbndsw)
      real(kind=8), intent(out) :: asma(nlay,nbndsw)         ! Aerosol asymmetry parameter
                                                      ! Dimensions: (nlay,nbndsw)

! Atmosphere/clouds - cldprop
      integer(kind=4), intent(out) :: inflag         ! flag for cloud property method
      integer(kind=4), intent(out) :: iceflag        ! flag for ice cloud properties
      integer(kind=4), intent(out) :: liqflag        ! flag for liquid cloud properties

      real(kind=8), intent(out) :: cldfrac(nlay)        ! layer cloud fraction
                                                      ! Dimensions: (nlay)
      real(kind=8), intent(out) :: tauc(nbndsw,nlay)         ! in-cloud optical depth (non-delta scaled)
                                                      ! Dimensions: (nbndsw,nlay)
      real(kind=8), intent(out) :: ssac(nbndsw,nlay)         ! in-cloud single scattering albedo (non-delta-scaled)
                                                      ! Dimensions: (nbndsw,nlay)
      real(kind=8), intent(out) :: asmc(nbndsw,nlay)         ! in-cloud asymmetry parameter (non-delta scaled)
                                                      ! Dimensions: (nbndsw,nlay)
      real(kind=8), intent(out) :: fsfc(nbndsw,nlay)         ! in-cloud forward scattering fraction (non-delta scaled)
                                                      ! Dimensions: (nbndsw,nlay)
      real(kind=8), intent(out) :: ciwp(nlay)           ! in-cloud ice water path
                                                      ! Dimensions: (nlay)
      real(kind=8), intent(out) :: clwp(nlay)           ! in-cloud liquid water path
                                                      ! Dimensions: (nlay)
      real(kind=8), intent(out) :: rel(nlay)            ! cloud liquid particle effective radius (microns)
                                                      ! Dimensions: (nlay)
      real(kind=8), intent(out) :: rei(nlay)            ! cloud ice particle effective size (microns)
                                                      ! Dimensions: (nlay)

! ----- Local -----
      real(kind=8), parameter :: amd = 28.9660    ! Effective molecular weight of dry air (g/mol)
      real(kind=8), parameter :: amw = 18.0160    ! Molecular weight of water vapor (g/mol)
!      real(kind=8), parameter :: amc = 44.0098   ! Molecular weight of carbon dioxide (g/mol)
!      real(kind=8), parameter :: amo = 47.9998   ! Molecular weight of ozone (g/mol)
!      real(kind=8), parameter :: amo2 = 31.9999  ! Molecular weight of oxygen (g/mol)
!      real(kind=8), parameter :: amch4 = 16.0430 ! Molecular weight of methane (g/mol)
!      real(kind=8), parameter :: amn2o = 44.0128 ! Molecular weight of nitrous oxide (g/mol)

! Set molecular weight ratios (for converting mmr to vmr)
!  e.g. h2ovmr = h2ommr * amdw)
      real(kind=8), parameter :: amdw = 1.607793  ! Molecular weight of dry air / water vapor
      real(kind=8), parameter :: amdc = 0.658114  ! Molecular weight of dry air / carbon dioxide
      real(kind=8), parameter :: amdo = 0.603428  ! Molecular weight of dry air / ozone
      real(kind=8), parameter :: amdm = 1.805423  ! Molecular weight of dry air / methane
      real(kind=8), parameter :: amdn = 0.658090  ! Molecular weight of dry air / nitrous oxide
      real(kind=8), parameter :: amdo2 = 0.905140 ! Molecular weight of dry air / oxygen

      real(kind=8), parameter :: sbc = 5.67e-08   ! Stefan-Boltzmann constant (W/m2K4)

      integer(kind=4) :: isp, l, ix, n, imol, ib       ! Loop indices
      real(kind=8) :: amm, summol                      ! 
      real(kind=8) :: adjflx                           ! flux adjustment for Earth/Sun distance
!      real(kind=8) :: earth_sun                        ! function for Earth/Sun distance adjustment

! Add one to nlayers here to include extra model layer at top of atmosphere
      nlayers = nlay

!  Initialize all molecular amounts to zero here, then pass input amounts
!  into RRTM array WKL below.

      wkl(:,:) = 0.0_8
      cldfrac(:) = 0.0_8
      tauc(:,:) = 0.0_8
      ssac(:,:) = 1.0_8
      asmc(:,:) = 0.0_8
      fsfc(:,:) = 0.0_8
      ciwp(:) = 0.0_8
      clwp(:) = 0.0_8
      rei(:) = 0.0_8
      rel(:) = 0.0_8
      taua(:,:) = 0.0_8
      ssaa(:,:) = 1.0_8
      asma(:,:) = 0.0_8
 
! Set flux adjustment for current Earth/Sun distance (two options).
! 1) Use Earth/Sun distance flux adjustment provided by GCM (input as adjes);
      adjflx = adjes
!
! 2) Calculate Earth/Sun distance from DYOFYR, the cumulative day of the year.
!    (Set adjflx to 1. to use constant Earth/Sun distance of 1 AU). 
      if (dyofyr .gt. 0) then
         adjflx = earth_sun(dyofyr)
      endif

! Set incoming solar flux adjustment to include adjustment for
! current Earth/Sun distance (ADJFLX) and scaling of default internal
! solar constant (rrsw_scon = 1368.22 Wm-2) by band (SOLVAR).  SOLVAR can be set 
! to a single scaling factor as needed, or to a different value in each 
! band, which may be necessary for paleoclimate simulations. 
! 
      do ib = jpb1,jpb2
!         solvar(ib) = 1._8
         solvar(ib) = scon / rrsw_scon
         adjflux(ib) = adjflx * solvar(ib)
      enddo

!  Set surface temperature.
      tbound = tsfc(iplon)

!  Install input GCM arrays into RRTMG_SW arrays for pressure, temperature,
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
         pdp(l) = pz(l-1) - pz(l)
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

! The following section can be used to set values for an additional layer (from
! the GCM top level to 1.e-4 mb) for improved calculation of TOA fluxes. 
! Temperature and molecular amounts in the extra model layer are set to 
! their values in the top GCM model layer, though these can be modified
! here if necessary. 
! If this feature is utilized, increase nlayers by one above, limit the two
! loops above to (nlayers-1), and set the top most (nlayers) layer values here. 

!      pavel(nlayers) = 0.5_8 * pz(nlayers-1)
!      tavel(nlayers) = tavel(nlayers-1)
!      pz(nlayers) = 1.e-4_8
!      tz(nlayers-1) = 0.5_8 * (tavel(nlayers)+tavel(nlayers-1))
!      tz(nlayers) = tz(nlayers-1)
!      pdp(nlayers) = pz(nlayers-1) - pz(nlayers)
!      wkl(1,nlayers) = wkl(1,nlayers-1)
!      wkl(2,nlayers) = wkl(2,nlayers-1)
!      wkl(3,nlayers) = wkl(3,nlayers-1)
!      wkl(4,nlayers) = wkl(4,nlayers-1)
!      wkl(6,nlayers) = wkl(6,nlayers-1)
!      wkl(7,nlayers) = wkl(7,nlayers-1)
!      amm = (1._8 - wkl(1,nlayers-1)) * amd + wkl(1,nlayers-1) * amw
!      coldry(nlayers) = (pz(nlayers-1)) * 1.e3_8 * avogad / &
!                        (1.e2_8 * grav * amm * (1._8 + wkl(1,nlayers-1)))

! At this point all molecular amounts in wkl are in volume mixing ratio; 
! convert to molec/cm2 based on coldry for use in rrtm.  

      do l = 1, nlayers
         do imol = 1, nmol
            wkl(imol,l) = coldry(l) * wkl(imol,l)
         enddo
      enddo

! Transfer aerosol optical properties to RRTM variables;
! modify to reverse layer indexing here if necessary.

      if (iaer .ge. 1) then 
         do l = 1, nlayers
            do ib = 1, nbndsw
               taua(l,ib) = tauaer(iplon,l,ib)
               ssaa(l,ib) = ssaaer(iplon,l,ib)
               asma(l,ib) = asmaer(iplon,l,ib)
            enddo
         enddo
      endif

! Transfer cloud fraction and cloud optical properties to RRTM variables;
! modify to reverse layer indexing here if necessary.

      if (icld .ge. 1) then 
         inflag = inflgsw
         iceflag = iceflgsw
         liqflag = liqflgsw

! Move incoming GCM cloud arrays to RRTMG cloud arrays.
! For GCM input, incoming reice is defined based on selected ice parameterization (inflglw)

         do l = 1, nlayers
            cldfrac(l) = cldfr(iplon,l)
            ciwp(l) = cicewp(iplon,l)
            clwp(l) = cliqwp(iplon,l)
            rei(l) = reice(iplon,l)
            rel(l) = reliq(iplon,l)
            do n = 1,nbndsw
               tauc(n,l) = taucld(n,iplon,l)
               ssac(n,l) = ssacld(n,iplon,l)
               asmc(n,l) = asmcld(n,iplon,l)
               fsfc(n,l) = fsfcld(n,iplon,l)
            enddo
         enddo

! If an extra layer is being used in RRTMG, set all cloud properties to zero in the extra layer.

!         cldfrac(nlayers) = 0.0_8
!         tauc(:,nlayers) = 0.0_8
!         ssac(:,nlayers) = 1.0_8
!         asmc(:,nlayers) = 0.0_8
!         fsfc(:,nlayers) = 0.0_8
!         ciwp(nlayers) = 0.0_8
!         clwp(nlayers) = 0.0_8
!         rei(nlayers) = 0.0_8
!         rel(nlayers) = 0.0_8
      
      endif

      end subroutine inatm_sw

      end module rrtmg_sw_rad


