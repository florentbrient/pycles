!     path:      $Source: /storm/rc1/cvsroot/rc/rrtmg_lw/src/rrtmg_lw_rtrn.f90,v $
!     author:    $Author: mike $
!     revision:  $Revision: 1.7 $
!     created:   $Date: 2009/11/12 20:52:25 $
!
      module rrtmg_lw_rtrn

!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2002-2009, Atmospheric & Environmental Research, Inc. (AER).  |
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
!  --------------------------------------------------------------------------

! --------- Modules ----------

      use parrrtm, only : mg, nbndlw, ngptlw
      use rrlw_con, only: fluxfac, heatfac
      use rrlw_wvn, only: delwave, ngs
      use rrlw_tbl, only: tblint, bpade, tau_tbl, exp_tbl, tfn_tbl
      use rrlw_vsn, only: hvrrtr, hnamrtr

      implicit none

      contains

!-----------------------------------------------------------------------------
      subroutine rtrn(nlayers, istart, iend, iout, pz, semiss, ncbands, &
                      cldfrac, taucloud, planklay, planklev, plankbnd, &
                      pwvcm, fracs, taut, & 
                      totuflux, totdflux, fnet, htr, &
                      totuclfl, totdclfl, fnetc, htrc, &
                      idrv, dplankbnd_dt, dtotuflux_dt, dtotuclfl_dt )
!-----------------------------------------------------------------------------
!
!  Original version:   E. J. Mlawer, et al. RRTM_V3.0
!  Revision for GCMs:  Michael J. Iacono; October, 2002
!  Revision for F90:  Michael J. Iacono; June, 2006
!  Revision for dFdT option: M. J. Iacono and E. J. Mlawer, November 2009
!
!  This program calculates the upward fluxes, downward fluxes, and
!  heating rates for an arbitrary clear or cloudy atmosphere.  The input
!  to this program is the atmospheric profile, all Planck function
!  information, and the cloud fraction by layer.  A variable diffusivity 
!  angle (SECDIFF) is used for the angle integration.  Bands 2-3 and 5-9 
!  use a value for SECDIFF that varies from 1.50 to 1.80 as a function of 
!  the column water vapor, and other bands use a value of 1.66.  The Gaussian 
!  weight appropriate to this angle (WTDIFF=0.5) is applied here.  Note that 
!  use of the emissivity angle for the flux integration can cause errors of 
!  1 to 4 W/m2 within cloudy layers.  
!  Clouds are treated with a random cloud overlap method.
!  This subroutine also provides the optional capability to calculate
!  the derivative of upward flux respect to surface temperature using
!  the pre-tabulated derivative of the Planck function with respect to 
!  temperature integrated over each spectral band.
!***************************************************************************

! ------- Declarations -------

! ----- Input -----
      integer(kind=4), intent(in) :: nlayers         ! total number of layers
      integer(kind=4), intent(in) :: istart          ! beginning band of calculation
      integer(kind=4), intent(in) :: iend            ! ending band of calculation
      integer(kind=4), intent(in) :: iout            ! output option flag

! Atmosphere
      real(kind=8), intent(in) :: pz(0:nlayers)             ! level (interface) pressures (hPa, mb)
                                                      !    Dimensions: (0:nlayers)
      real(kind=8), intent(in) :: pwvcm              ! precipitable water vapor (cm)
      real(kind=8), intent(in) :: semiss(nbndlw)          ! lw surface emissivity
                                                      !    Dimensions: (nbndlw)
      real(kind=8), intent(in) :: planklay(nlayers,nbndlw)      ! 
                                                      !    Dimensions: (nlayers,nbndlw)
      real(kind=8), intent(in) :: planklev(0:nlayers,nbndlw)     ! 
                                                      !    Dimensions: (0:nlayers,nbndlw)
      real(kind=8), intent(in) :: plankbnd(nbndlw)        ! 
                                                      !    Dimensions: (nbndlw)
      real(kind=8), intent(in) :: fracs(nlayers,ngptlw)         ! 
                                                      !    Dimensions: (nlayers,ngptw)
      real(kind=8), intent(in) :: taut(nlayers,ngptlw)          ! gaseous + aerosol optical depths
                                                      !    Dimensions: (nlayers,ngptlw)

! Clouds
      integer(kind=4), intent(in) :: ncbands         ! number of cloud spectral bands
      real(kind=8), intent(in) :: cldfrac(nlayers)         ! layer cloud fraction
                                                      !    Dimensions: (nlayers)
      real(kind=8), intent(in) :: taucloud(nlayers,nbndlw)      ! layer cloud optical depth
                                                      !    Dimensions: (nlayers,nbndlw)
      integer(kind=4), intent(in) :: idrv            ! flag for calculation of dF/dt from 
                                                      ! Planck derivative [0=off, 1=on]
      real(kind=8), intent(in) :: dplankbnd_dt(nbndlw)    ! derivative of Planck function wrt temp
                                                      !    Dimensions: (nbndlw)

! ----- Output -----
      real(kind=8), intent(out) :: totuflux(0:nlayers)      ! upward longwave flux (w/m2)
                                                      !    Dimensions: (0:nlayers)
      real(kind=8), intent(out) :: totdflux(0:nlayers)      ! downward longwave flux (w/m2)
                                                      !    Dimensions: (0:nlayers)
      real(kind=8), intent(out) :: fnet(0:nlayers)          ! net longwave flux (w/m2)
                                                      !    Dimensions: (0:nlayers)
      real(kind=8), intent(out) :: htr(0:nlayers)           ! longwave heating rate (k/day)
                                                      !    Dimensions: (0:nlayers)
      real(kind=8), intent(out) :: totuclfl(0:nlayers)      ! clear sky upward longwave flux (w/m2)
                                                      !    Dimensions: (0:nlayers)
      real(kind=8), intent(out) :: totdclfl(0:nlayers)      ! clear sky downward longwave flux (w/m2)
                                                      !    Dimensions: (0:nlayers)
      real(kind=8), intent(out) :: fnetc(0:nlayers)         ! clear sky net longwave flux (w/m2)
                                                      !    Dimensions: (0:nlayers)
      real(kind=8), intent(out) :: htrc(0:nlayers)          ! clear sky longwave heating rate (k/day)
                                                      !    Dimensions: (0:nlayers)
      real(kind=8), intent(out) :: dtotuflux_dt(0:nlayers)  ! change in upward longwave flux (w/m2/k)
                                                      ! with respect to surface temperature
                                                      !    Dimensions: (0:nlayers)
      real(kind=8), intent(out) :: dtotuclfl_dt(0:nlayers)  ! change in upward longwave flux (w/m2/k)
                                                      ! with respect to surface temperature
                                                      !    Dimensions: (0:nlayers)

! ----- Local -----
! Declarations for radiative transfer
      real(kind=8) :: abscld(nlayers,nbndlw)
      real(kind=8) :: atot(nlayers)
      real(kind=8) :: atrans(nlayers)
      real(kind=8) :: bbugas(nlayers)
      real(kind=8) :: bbutot(nlayers)
      real(kind=8) :: clrurad(0:nlayers)
      real(kind=8) :: clrdrad(0:nlayers)
      real(kind=8) :: efclfrac(nlayers,nbndlw)
      real(kind=8) :: uflux(0:nlayers)
      real(kind=8) :: dflux(0:nlayers)
      real(kind=8) :: urad(0:nlayers)
      real(kind=8) :: drad(0:nlayers)
      real(kind=8) :: uclfl(0:nlayers)
      real(kind=8) :: dclfl(0:nlayers)
      real(kind=8) :: odcld(nlayers,nbndlw)


      real(kind=8) :: secdiff(nbndlw)                 ! secant of diffusivity angle
      real(kind=8) :: a0(nbndlw),a1(nbndlw),a2(nbndlw)! diffusivity angle adjustment coefficients
      real(kind=8) :: wtdiff, rec_6
      real(kind=8) :: transcld, radld, radclrd, plfrac, blay, dplankup, dplankdn
      real(kind=8) :: odepth, odtot, odepth_rec, odtot_rec, gassrc
      real(kind=8) :: tblind, tfactot, bbd, bbdtot, tfacgas, transc, tausfac
      real(kind=8) :: rad0, reflect, radlu, radclru

      real(kind=8) :: duflux_dt(0:nlayers)
      real(kind=8) :: duclfl_dt(0:nlayers)
      real(kind=8) :: d_urad_dt(0:nlayers)
      real(kind=8) :: d_clrurad_dt(0:nlayers)
      real(kind=8) :: d_rad0_dt, d_radlu_dt, d_radclru_dt

      integer(kind=4) :: icldlyr(nlayers)             ! flag for cloud in layer
      integer(kind=4) :: ibnd, ib, iband, lay, lev, l ! loop indices
      integer(kind=4) :: igc                          ! g-point interval counter
      integer(kind=4) :: iclddn                       ! flag for cloud in down path
      integer(kind=4) :: ittot, itgas, itr            ! lookup table indices
      integer(kind=4) :: ipat(16,0:2)


! ------- Definitions -------
! input
!    nlayers                      ! number of model layers
!    ngptlw                       ! total number of g-point subintervals
!    nbndlw                       ! number of longwave spectral bands
!    ncbands                      ! number of spectral bands for clouds
!    secdiff                      ! diffusivity angle
!    wtdiff                       ! weight for radiance to flux conversion
!    pavel                        ! layer pressures (mb)
!    pz                           ! level (interface) pressures (mb)
!    tavel                        ! layer temperatures (k)
!    tz                           ! level (interface) temperatures(mb)
!    tbound                       ! surface temperature (k)
!    cldfrac                      ! layer cloud fraction
!    taucloud                     ! layer cloud optical depth
!    itr                          ! integer look-up table index
!    icldlyr                      ! flag for cloudy layers
!    iclddn                       ! flag for cloud in column at any layer
!    semiss                       ! surface emissivities for each band
!    reflect                      ! surface reflectance
!    bpade                        ! 1/(pade constant)
!    tau_tbl                      ! clear sky optical depth look-up table
!    exp_tbl                      ! exponential look-up table for transmittance
!    tfn_tbl                      ! tau transition function look-up table

! local
!    atrans                       ! gaseous absorptivity
!    abscld                       ! cloud absorptivity
!    atot                         ! combined gaseous and cloud absorptivity
!    odclr                        ! clear sky (gaseous) optical depth
!    odcld                        ! cloud optical depth
!    odtot                        ! optical depth of gas and cloud
!    tfacgas                      ! gas-only pade factor, used for planck fn
!    tfactot                      ! gas and cloud pade factor, used for planck fn
!    bbdgas                       ! gas-only planck function for downward rt
!    bbugas                       ! gas-only planck function for upward rt
!    bbdtot                       ! gas and cloud planck function for downward rt
!    bbutot                       ! gas and cloud planck function for upward calc.
!    gassrc                       ! source radiance due to gas only
!    efclfrac                     ! effective cloud fraction
!    radlu                        ! spectrally summed upward radiance 
!    radclru                      ! spectrally summed clear sky upward radiance 
!    urad                         ! upward radiance by layer
!    clrurad                      ! clear sky upward radiance by layer
!    radld                        ! spectrally summed downward radiance 
!    radclrd                      ! spectrally summed clear sky downward radiance 
!    drad                         ! downward radiance by layer
!    clrdrad                      ! clear sky downward radiance by layer
!    d_radlu_dt                   ! spectrally summed upward radiance 
!    d_radclru_dt                 ! spectrally summed clear sky upward radiance 
!    d_urad_dt                    ! upward radiance by layer
!    d_clrurad_dt                 ! clear sky upward radiance by layer

! output
!    totuflux                     ! upward longwave flux (w/m2)
!    totdflux                     ! downward longwave flux (w/m2)
!    fnet                         ! net longwave flux (w/m2)
!    htr                          ! longwave heating rate (k/day)
!    totuclfl                     ! clear sky upward longwave flux (w/m2)
!    totdclfl                     ! clear sky downward longwave flux (w/m2)
!    fnetc                        ! clear sky net longwave flux (w/m2)
!    htrc                         ! clear sky longwave heating rate (k/day)
!    dtotuflux_dt                 ! change in upward longwave flux (w/m2/k)
!                                 ! with respect to surface temperature
!    dtotuclfl_dt                 ! change in clear sky upward longwave flux (w/m2/k)
!                                 ! with respect to surface temperature

! These arrays indicate the spectral 'region' (used in the 
! calculation of ice cloud optical depths) corresponding
! to each spectral band.  See cldprop.f for more details.
      data ipat /1,1,1,1,1,1,1,1,1, 1, 1, 1, 1, 1, 1, 1, &
                 1,2,3,3,3,4,4,4,5, 5, 5, 5, 5, 5, 5, 5, &
                 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16/

! This secant and weight corresponds to the standard diffusivity 
! angle.  This initial value is redefined below for some bands.
      data wtdiff /0.5_8/
      data rec_6 /0.166667_8/

! Reset diffusivity angle for Bands 2-3 and 5-9 to vary (between 1.50
! and 1.80) as a function of total column water vapor.  The function
! has been defined to minimize flux and cooling rate errors in these bands
! over a wide range of precipitable water values.
      data a0 / 1.66_8,  1.55_8,  1.58_8,  1.66_8, &
                1.54_8, 1.454_8,  1.89_8,  1.33_8, &
               1.668_8,  1.66_8,  1.66_8,  1.66_8, &
                1.66_8,  1.66_8,  1.66_8,  1.66_8 /
      data a1 / 0.00_8,  0.25_8,  0.22_8,  0.00_8, &
                0.13_8, 0.446_8, -0.10_8,  0.40_8, &
              -0.006_8,  0.00_8,  0.00_8,  0.00_8, &
                0.00_8,  0.00_8,  0.00_8,  0.00_8 /
      data a2 / 0.00_8, -12.0_8, -11.7_8,  0.00_8, &
               -0.72_8,-0.243_8,  0.19_8,-0.062_8, &
               0.414_8,  0.00_8,  0.00_8,  0.00_8, &
                0.00_8,  0.00_8,  0.00_8,  0.00_8 /

      hvrrtr = '$Revision: 1.7 $'

      do ibnd = 1,nbndlw
         if (ibnd.eq.1 .or. ibnd.eq.4 .or. ibnd.ge.10) then
           secdiff(ibnd) = 1.66_8
         else
           secdiff(ibnd) = a0(ibnd) + a1(ibnd)*exp(a2(ibnd)*pwvcm)
           if (secdiff(ibnd) .gt. 1.80_8) secdiff(ibnd) = 1.80_8
           if (secdiff(ibnd) .lt. 1.50_8) secdiff(ibnd) = 1.50_8
         endif
      enddo

      urad(0) = 0.0_8
      drad(0) = 0.0_8
      totuflux(0) = 0.0_8
      totdflux(0) = 0.0_8
      clrurad(0) = 0.0_8
      clrdrad(0) = 0.0_8
      totuclfl(0) = 0.0_8
      totdclfl(0) = 0.0_8
      if (idrv .eq. 1) then
         d_urad_dt(0) = 0.0_8
         d_clrurad_dt(0) = 0.0_8
         dtotuflux_dt(0) = 0.0_8
         dtotuclfl_dt(0) = 0.0_8
      endif

      do lay = 1, nlayers
         urad(lay) = 0.0_8
         drad(lay) = 0.0_8
         totuflux(lay) = 0.0_8
         totdflux(lay) = 0.0_8
         clrurad(lay) = 0.0_8
         clrdrad(lay) = 0.0_8
         totuclfl(lay) = 0.0_8
         totdclfl(lay) = 0.0_8
         if (idrv .eq. 1) then
            d_urad_dt(lay) = 0.0_8
            d_clrurad_dt(lay) = 0.0_8
            dtotuflux_dt(lay) = 0.0_8
            dtotuclfl_dt(lay) = 0.0_8
         endif

         do ib = 1, ncbands
            if (cldfrac(lay) .ge. 1.e-6_8) then
               odcld(lay,ib) = secdiff(ib) * taucloud(lay,ib)
               transcld = exp(-odcld(lay,ib))
               abscld(lay,ib) = 1. - transcld
               efclfrac(lay,ib) = abscld(lay,ib) * cldfrac(lay)
               icldlyr(lay) = 1
            else
               odcld(lay,ib) = 0.0_8
               abscld(lay,ib) = 0.0_8
               efclfrac(lay,ib) = 0.0_8
               icldlyr(lay) = 0
            endif
         enddo
      enddo

      igc = 1
! Loop over frequency bands.
      do iband = istart, iend

! Reinitialize g-point counter for each band if output for each band is requested.
         if (iout.gt.0.and.iband.ge.2) igc = ngs(iband-1)+1
         if (ncbands .eq. 1) then
            ib = ipat(iband,0)
         elseif (ncbands .eq.  5) then
            ib = ipat(iband,1)
         elseif (ncbands .eq. 16) then
            ib = ipat(iband,2)
         endif

! Loop over g-channels.
 1000    continue

! Radiative transfer starts here.
         radld = 0._8
         radclrd = 0._8
         iclddn = 0

! Downward radiative transfer loop.  

         do lev = nlayers, 1, -1
               plfrac = fracs(lev,igc)
               blay = planklay(lev,iband)
               dplankup = planklev(lev,iband) - blay
               dplankdn = planklev(lev-1,iband) - blay
               odepth = secdiff(iband) * taut(lev,igc)
               if (odepth .lt. 0.0_8) odepth = 0.0_8
! Cloudy layer
               if (icldlyr(lev).eq.1) then
                  iclddn = 1
                  odtot = odepth + odcld(lev,ib)
                  if (odtot .lt. 0.06_8) then
                     atrans(lev) = odepth - 0.5_8*odepth*odepth
                     odepth_rec = rec_6*odepth
                     gassrc = plfrac*(blay+dplankdn*odepth_rec)*atrans(lev)

                     atot(lev) =  odtot - 0.5_8*odtot*odtot
                     odtot_rec = rec_6*odtot
                     bbdtot =  plfrac * (blay+dplankdn*odtot_rec)
                     bbd = plfrac*(blay+dplankdn*odepth_rec)
                     radld = radld - radld * (atrans(lev) + &
                         efclfrac(lev,ib) * (1. - atrans(lev))) + &
                         gassrc + cldfrac(lev) * &
                         (bbdtot * atot(lev) - gassrc)
                     drad(lev-1) = drad(lev-1) + radld
                  
                     bbugas(lev) =  plfrac * (blay+dplankup*odepth_rec)
                     bbutot(lev) =  plfrac * (blay+dplankup*odtot_rec)

                  elseif (odepth .le. 0.06_8) then
                     atrans(lev) = odepth - 0.5_8*odepth*odepth
                     odepth_rec = rec_6*odepth
                     gassrc = plfrac*(blay+dplankdn*odepth_rec)*atrans(lev)

                     odtot = odepth + odcld(lev,ib)
                     tblind = odtot/(bpade+odtot)
                     ittot = tblint*tblind + 0.5_8
                     tfactot = tfn_tbl(ittot)
                     bbdtot = plfrac * (blay + tfactot*dplankdn)
                     bbd = plfrac*(blay+dplankdn*odepth_rec)
                     atot(lev) = 1._8 - exp_tbl(ittot)

                     radld = radld - radld * (atrans(lev) + &
                         efclfrac(lev,ib) * (1._8 - atrans(lev))) + &
                         gassrc + cldfrac(lev) * &
                         (bbdtot * atot(lev) - gassrc)
                     drad(lev-1) = drad(lev-1) + radld

                     bbugas(lev) = plfrac * (blay + dplankup*odepth_rec)
                     bbutot(lev) = plfrac * (blay + tfactot * dplankup)

                  else

                     tblind = odepth/(bpade+odepth)
                     itgas = tblint*tblind+0.5_8
                     odepth = tau_tbl(itgas)
                     atrans(lev) = 1._8 - exp_tbl(itgas)
                     tfacgas = tfn_tbl(itgas)
                     gassrc = atrans(lev) * plfrac * (blay + tfacgas*dplankdn)

                     odtot = odepth + odcld(lev,ib)
                     tblind = odtot/(bpade+odtot)
                     ittot = tblint*tblind + 0.5_8
                     tfactot = tfn_tbl(ittot)
                     bbdtot = plfrac * (blay + tfactot*dplankdn)
                     bbd = plfrac*(blay+tfacgas*dplankdn)
                     atot(lev) = 1._8 - exp_tbl(ittot)

                  radld = radld - radld * (atrans(lev) + &
                     efclfrac(lev,ib) * (1._8 - atrans(lev))) + &
                     gassrc + cldfrac(lev) * &
                     (bbdtot * atot(lev) - gassrc)
                  drad(lev-1) = drad(lev-1) + radld
                  bbugas(lev) = plfrac * (blay + tfacgas * dplankup)
                  bbutot(lev) = plfrac * (blay + tfactot * dplankup)
                  endif
! Clear layer
               else
                  if (odepth .le. 0.06_8) then
                     atrans(lev) = odepth-0.5_8*odepth*odepth
                     odepth = rec_6*odepth
                     bbd = plfrac*(blay+dplankdn*odepth)
                     bbugas(lev) = plfrac*(blay+dplankup*odepth)
                  else
                     tblind = odepth/(bpade+odepth)
                     itr = tblint*tblind+0.5_8
                     transc = exp_tbl(itr)
                     atrans(lev) = 1._8-transc
                     tausfac = tfn_tbl(itr)
                     bbd = plfrac*(blay+tausfac*dplankdn)
                     bbugas(lev) = plfrac * (blay + tausfac * dplankup)
                  endif   
                  radld = radld + (bbd-radld)*atrans(lev)
                  drad(lev-1) = drad(lev-1) + radld
               endif
!  Set clear sky stream to total sky stream as long as layers
!  remain clear.  Streams diverge when a cloud is reached (iclddn=1),
!  and clear sky stream must be computed separately from that point.
                  if (iclddn.eq.1) then
                     radclrd = radclrd + (bbd-radclrd) * atrans(lev) 
                     clrdrad(lev-1) = clrdrad(lev-1) + radclrd
                  else
                     radclrd = radld
                     clrdrad(lev-1) = drad(lev-1)
                  endif
            enddo

! Spectral emissivity & reflectance
!  Include the contribution of spectrally varying longwave emissivity
!  and reflection from the surface to the upward radiative transfer.
!  Note: Spectral and Lambertian reflection are identical for the
!  diffusivity angle flux integration used here.
!  Note: The emissivity is applied to plankbnd and dplankbnd_dt when 
!  they are defined in subroutine setcoef. 

         rad0 = fracs(1,igc) * plankbnd(iband)
         if (idrv .eq. 1) then
            d_rad0_dt = fracs(1,igc) * dplankbnd_dt(iband)
         endif

!  Add in specular reflection of surface downward radiance.
         reflect = 1._8 - semiss(iband)
         radlu = rad0 + reflect * radld
         radclru = rad0 + reflect * radclrd


! Upward radiative transfer loop.
         urad(0) = urad(0) + radlu
         clrurad(0) = clrurad(0) + radclru
         if (idrv .eq. 1) then
            d_radlu_dt = d_rad0_dt
            d_urad_dt(0) = d_urad_dt(0) + d_radlu_dt
            d_radclru_dt = d_rad0_dt
            d_clrurad_dt(0) = d_clrurad_dt(0) + d_radclru_dt
         endif

         do lev = 1, nlayers
! Cloudy layer
            if (icldlyr(lev) .eq. 1) then
               gassrc = bbugas(lev) * atrans(lev)
               radlu = radlu - radlu * (atrans(lev) + &
                   efclfrac(lev,ib) * (1._8 - atrans(lev))) + &
                   gassrc + cldfrac(lev) * &
                   (bbutot(lev) * atot(lev) - gassrc)
               urad(lev) = urad(lev) + radlu
               if (idrv .eq. 1) then
                  d_radlu_dt = d_radlu_dt * cldfrac(lev) * (1.0_8 - atot(lev)) + &
                         d_radlu_dt * (1.0_8 - cldfrac(lev)) * (1.0_8 - atrans(lev))
                  d_urad_dt(lev) = d_urad_dt(lev) + d_radlu_dt
               endif
! Clear layer
            else
               radlu = radlu + (bbugas(lev)-radlu)*atrans(lev)
               urad(lev) = urad(lev) + radlu
               if (idrv .eq. 1) then
                  d_radlu_dt = d_radlu_dt * (1.0_8 - atrans(lev))
                  d_urad_dt(lev) = d_urad_dt(lev) + d_radlu_dt
               endif
            endif
!  Set clear sky stream to total sky stream as long as all layers
!  are clear (iclddn=0).  Streams must be calculated separately at 
!  all layers when a cloud is present (iclddn=1), because surface 
!  reflectance is different for each stream.
               if (iclddn.eq.1) then
                  radclru = radclru + (bbugas(lev)-radclru)*atrans(lev) 
                  clrurad(lev) = clrurad(lev) + radclru
               else
                  radclru = radlu
                  clrurad(lev) = urad(lev)
               endif
               if (idrv .eq. 1) then
                  if (iclddn.eq.1) then
                     d_radclru_dt = d_radclru_dt * (1.0_8 - atrans(lev))
                     d_clrurad_dt(lev) = d_clrurad_dt(lev) + d_radclru_dt
                  else
                     d_radclru_dt = d_radlu_dt
                     d_clrurad_dt(lev) = d_urad_dt(lev)
                  endif
               endif
         enddo

! Increment g-point counter
         igc = igc + 1
! Return to continue radiative transfer for all g-channels in present band
         if (igc .le. ngs(iband)) go to 1000

! Process longwave output from band for total and clear streams.
! Calculate upward, downward, and net flux.
         do lev = nlayers, 0, -1
            uflux(lev) = urad(lev)*wtdiff
            dflux(lev) = drad(lev)*wtdiff
            urad(lev) = 0.0_8
            drad(lev) = 0.0_8
            totuflux(lev) = totuflux(lev) + uflux(lev) * delwave(iband)
            totdflux(lev) = totdflux(lev) + dflux(lev) * delwave(iband)
            uclfl(lev) = clrurad(lev)*wtdiff
            dclfl(lev) = clrdrad(lev)*wtdiff
            clrurad(lev) = 0.0_8
            clrdrad(lev) = 0.0_8
            totuclfl(lev) = totuclfl(lev) + uclfl(lev) * delwave(iband)
            totdclfl(lev) = totdclfl(lev) + dclfl(lev) * delwave(iband)
         enddo

! Calculate total change in upward flux wrt surface temperature
         if (idrv .eq. 1) then
            do lev = nlayers, 0, -1
               duflux_dt(lev) = d_urad_dt(lev) * wtdiff
               d_urad_dt(lev) = 0.0_8
               dtotuflux_dt(lev) = dtotuflux_dt(lev) + duflux_dt(lev) * delwave(iband) * fluxfac
               duclfl_dt(lev) = d_clrurad_dt(lev) * wtdiff
               d_clrurad_dt(lev) = 0.0_8
               dtotuclfl_dt(lev) = dtotuclfl_dt(lev) + duclfl_dt(lev) * delwave(iband) * fluxfac
            enddo
         endif

! End spectral band loop
      enddo

! Calculate fluxes at surface
      totuflux(0) = totuflux(0) * fluxfac
      totdflux(0) = totdflux(0) * fluxfac
      fnet(0) = totuflux(0) - totdflux(0)
      totuclfl(0) = totuclfl(0) * fluxfac
      totdclfl(0) = totdclfl(0) * fluxfac
      fnetc(0) = totuclfl(0) - totdclfl(0)

! Calculate fluxes at model levels
      do lev = 1, nlayers
         totuflux(lev) = totuflux(lev) * fluxfac
         totdflux(lev) = totdflux(lev) * fluxfac
         fnet(lev) = totuflux(lev) - totdflux(lev)
         totuclfl(lev) = totuclfl(lev) * fluxfac
         totdclfl(lev) = totdclfl(lev) * fluxfac
         fnetc(lev) = totuclfl(lev) - totdclfl(lev)
         l = lev - 1

! Calculate heating rates at model layers
         htr(l)=heatfac*(fnet(l)-fnet(lev))/(pz(l)-pz(lev)) 
         htrc(l)=heatfac*(fnetc(l)-fnetc(lev))/(pz(l)-pz(lev)) 
      enddo

! Set heating rate to zero in top layer
      htr(nlayers) = 0.0_8
      htrc(nlayers) = 0.0_8

      end subroutine rtrn

      end module rrtmg_lw_rtrn

