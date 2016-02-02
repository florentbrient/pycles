!     path:      $Source$
!     author:    $Author: mike $
!     revision:  $Revision: 11661 $
!     created:   $Date: 2009-05-22 18:22:22 -0400 (Fri, 22 May 2009) $

      module rrtmg_sw_spcvrt

!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2002-2009, Atmospheric & Environmental Research, Inc. (AER).  |
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
!  --------------------------------------------------------------------------

! ------- Modules -------

      
      use parrrsw, only : nbndsw, ngptsw, mxmol, jpband
      use rrsw_tbl, only : tblint, bpade, od_lo, exp_tbl
      use rrsw_vsn, only : hvrspv, hnamspv
      use rrsw_wvn, only : ngc, ngs
      use rrtmg_sw_reftra, only: reftra_sw
      use rrtmg_sw_taumol, only: taumol_sw
      use rrtmg_sw_vrtqdr, only: vrtqdr_sw

      implicit none

      contains

! ---------------------------------------------------------------------------
      subroutine spcvrt_sw &
            (nlayers, istart, iend, icpr, idelm, iout, &
             pavel, tavel, pz, tz, tbound, palbd, palbp, &
             pclfr, ptauc, pasyc, pomgc, ptaucorig, &
             ptaua, pasya, pomga, prmu0, coldry, wkl, adjflux, &
             laytrop, layswtch, laylow, jp, jt, jt1, &
             co2mult, colch4, colco2, colh2o, colmol, coln2o, colo2, colo3, &
             fac00, fac01, fac10, fac11, &
             selffac, selffrac, indself, forfac, forfrac, indfor, &
             pbbfd, pbbfu, pbbcd, pbbcu, puvfd, puvcd, pnifd, pnicd, &
             pbbfddir, pbbcddir, puvfddir, puvcddir, pnifddir, pnicddir)
! ---------------------------------------------------------------------------
!
! Purpose: Contains spectral loop to compute the shortwave radiative fluxes, 
!          using the two-stream method of H. Barker. 
!
! Interface:  *spcvrt_sw* is called from *rrtmg_sw.F90* or rrtmg_sw.1col.F90*
!
! Method:
!    Adapted from two-stream model of H. Barker;
!    Two-stream model options (selected with kmodts in rrtmg_sw_reftra.F90):
!        1: Eddington, 2: PIFM, Zdunkowski et al., 3: discret ordinates
!
! Modifications:
!
! Original: H. Barker
! Revision: Merge with RRTMG_SW: J.-J.Morcrette, ECMWF, Feb 2003
! Revision: Add adjustment for Earth/Sun distance : MJIacono, AER, Oct 2003
! Revision: Bug fix for use of PALBP and PALBD: MJIacono, AER, Nov 2003
! Revision: Bug fix to apply delta scaling to clear sky: AER, Dec 2004
! Revision: Code modified so that delta scaling is not done in cloudy profiles
!           if routine cldprop is used; delta scaling can be applied by swithcing
!           code below if cldprop is not used to get cloud properties. 
!           AER, Jan 2005
! Revision: Uniform formatting for RRTMG: MJIacono, AER, Jul 2006 
! Revision: Use exponential lookup table for transmittance: MJIacono, AER, 
!           Aug 2007 
!
! ------------------------------------------------------------------

! ------- Declarations ------

! -------- Input -------

      integer(kind=4), intent(in) :: nlayers
      integer(kind=4), intent(in) :: istart
      integer(kind=4), intent(in) :: iend
      integer(kind=4), intent(in) :: icpr
      integer(kind=4), intent(in) :: idelm   ! delta-m scaling flag
                                              ! [0 = direct and diffuse fluxes are unscaled]
                                              ! [1 = direct and diffuse fluxes are scaled]
      integer(kind=4), intent(in) :: iout
      integer(kind=4), intent(in) :: laytrop
      integer(kind=4), intent(in) :: layswtch
      integer(kind=4), intent(in) :: laylow

      integer(kind=4), intent(in) :: indfor(nlayers)
                                                               !   Dimensions: (nlayers)
      integer(kind=4), intent(in) :: indself(nlayers)
                                                               !   Dimensions: (nlayers)
      integer(kind=4), intent(in) :: jp(nlayers)
                                                               !   Dimensions: (nlayers)
      integer(kind=4), intent(in) :: jt(nlayers)
                                                               !   Dimensions: (nlayers)
      integer(kind=4), intent(in) :: jt1(nlayers)
                                                               !   Dimensions: (nlayers)

      real(kind=8), intent(in) :: pavel(nlayers)                    ! layer pressure (hPa, mb) 
                                                               !   Dimensions: (nlayers)
      real(kind=8), intent(in) :: tavel(nlayers)                    ! layer temperature (K)
                                                               !   Dimensions: (nlayers)
      real(kind=8), intent(in) :: pz(0:nlayers)                      ! level (interface) pressure (hPa, mb)
                                                               !   Dimensions: (0:nlayers)
      real(kind=8), intent(in) :: tz(0:nlayers)                      ! level temperatures (hPa, mb)
                                                               !   Dimensions: (0:nlayers)
      real(kind=8), intent(in) :: tbound                      ! surface temperature (K)
      real(kind=8), intent(in) :: wkl(mxmol,nlayers)                    ! molecular amounts (mol/cm2) 
                                                               !   Dimensions: (mxmol,nlayers)
      real(kind=8), intent(in) :: coldry(nlayers)                   ! dry air column density (mol/cm2)
                                                               !   Dimensions: (nlayers)
      real(kind=8), intent(in) :: colmol(nlayers)
                                                               !   Dimensions: (nlayers)
      real(kind=8), intent(in) :: adjflux(jpband)                  ! Earth/Sun distance adjustment
                                                               !   Dimensions: (jpband)

      real(kind=8), intent(in) :: palbd(nbndsw)                    ! surface albedo (diffuse)
                                                               !   Dimensions: (nbndsw)
      real(kind=8), intent(in) :: palbp(nbndsw)                    ! surface albedo (direct)
                                                               !   Dimensions: (nbndsw)
      real(kind=8), intent(in) :: prmu0                       ! cosine of solar zenith angle
      real(kind=8), intent(in) :: pclfr(nlayers)                    ! cloud fraction
                                                               !   Dimensions: (nlayers)
      real(kind=8), intent(in) :: ptauc(nlayers,nbndsw)                  ! cloud optical depth
                                                               !   Dimensions: (nlayers,nbndsw)
      real(kind=8), intent(in) :: pasyc(nlayers,nbndsw)                  ! cloud asymmetry parameter
                                                               !   Dimensions: (nlayers,nbndsw)
      real(kind=8), intent(in) :: pomgc(nlayers,nbndsw)                  ! cloud single scattering albedo
                                                               !   Dimensions: (nlayers,nbndsw)
      real(kind=8), intent(in) :: ptaucorig(nlayers,nbndsw)              ! cloud optical depth, non-delta scaled
                                                               !   Dimensions: (nlayers,nbndsw)
      real(kind=8), intent(in) :: ptaua(nlayers,nbndsw)                  ! aerosol optical depth
                                                               !   Dimensions: (nlayers,nbndsw)
      real(kind=8), intent(in) :: pasya(nlayers,nbndsw)                  ! aerosol asymmetry parameter
                                                               !   Dimensions: (nlayers,nbndsw)
      real(kind=8), intent(in) :: pomga(nlayers,nbndsw)                  ! aerosol single scattering albedo
                                                               !   Dimensions: (nlayers,nbndsw)

      real(kind=8), intent(in) :: colh2o(nlayers)
                                                               !   Dimensions: (nlayers)
      real(kind=8), intent(in) :: colco2(nlayers)
                                                               !   Dimensions: (nlayers)
      real(kind=8), intent(in) :: colch4(nlayers)
                                                               !   Dimensions: (nlayers)
      real(kind=8), intent(in) :: co2mult(nlayers)
                                                               !   Dimensions: (nlayers)
      real(kind=8), intent(in) :: colo3(nlayers)
                                                               !   Dimensions: (nlayers)
      real(kind=8), intent(in) :: colo2(nlayers)
                                                               !   Dimensions: (nlayers)
      real(kind=8), intent(in) :: coln2o(nlayers)
                                                               !   Dimensions: (nlayers)

      real(kind=8), intent(in) :: forfac(nlayers)
                                                               !   Dimensions: (nlayers)
      real(kind=8), intent(in) :: forfrac(nlayers)
                                                               !   Dimensions: (nlayers)
      real(kind=8), intent(in) :: selffac(nlayers)
                                                               !   Dimensions: (nlayers)
      real(kind=8), intent(in) :: selffrac(nlayers)
                                                               !   Dimensions: (nlayers)
      real(kind=8), intent(in) :: fac00(nlayers)
                                                               !   Dimensions: (nlayers)
      real(kind=8), intent(in) :: fac01(nlayers)
                                                               !   Dimensions: (nlayers)
      real(kind=8), intent(in) :: fac10(nlayers)
                                                               !   Dimensions: (nlayers)
      real(kind=8), intent(in) :: fac11(nlayers)
                                                               !   Dimensions: (nlayers)

! ------- Output -------
                                                               !   All Dimensions: (nlayers+1)
      real(kind=8), intent(out) :: pbbcd(nlayers+1)
      real(kind=8), intent(out) :: pbbcu(nlayers+1)
      real(kind=8), intent(out) :: pbbfd(nlayers+1)
      real(kind=8), intent(out) :: pbbfu(nlayers+1)
      real(kind=8), intent(out) :: pbbfddir(nlayers+1)
      real(kind=8), intent(out) :: pbbcddir(nlayers+1)

      real(kind=8), intent(out) :: puvcd(nlayers+1)
      real(kind=8), intent(out) :: puvfd(nlayers+1)
      real(kind=8), intent(out) :: puvcddir(nlayers+1)
      real(kind=8), intent(out) :: puvfddir(nlayers+1)

      real(kind=8), intent(out) :: pnicd(nlayers+1)
      real(kind=8), intent(out) :: pnifd(nlayers+1)
      real(kind=8), intent(out) :: pnicddir(nlayers+1)
      real(kind=8), intent(out) :: pnifddir(nlayers+1)

! Output - inactive                                            !   All Dimensions: (nlayers+1)
!      real(kind=8), intent(out) :: puvcu(nlayers+1)
!      real(kind=8), intent(out) :: puvfu(nlayers+1)
!      real(kind=8), intent(out) :: pnicu(nlayers+1)
!      real(kind=8), intent(out) :: pnifu(nlayers+1)
!      real(kind=8), intent(out) :: pvscd(nlayers+1)
!      real(kind=8), intent(out) :: pvscu(nlayers+1)
!      real(kind=8), intent(out) :: pvsfd(nlayers+1)
!      real(kind=8), intent(out) :: pvsfu(nlayers+1)


! ------- Local -------

      logical :: lrtchkclr(nlayers),lrtchkcld(nlayers)

      integer(kind=4)  :: klev
      integer(kind=4) :: ib1, ib2, ibm, igt, ikl, ikp, ikx
      integer(kind=4) :: iw, jb, jg, jl, jk
!      integer(kind=4), parameter :: nuv = ?? 
!      integer(kind=4), parameter :: nvs = ?? 
      integer(kind=4) :: itind

      real(kind=8) :: tblind, ze1
      real(kind=8) :: zclear, zcloud
      real(kind=8) :: zdbt(nlayers+1), zdbt_nodel(nlayers+1)
      real(kind=8) :: zgc(nlayers), zgcc(nlayers), zgco(nlayers)
      real(kind=8) :: zomc(nlayers), zomcc(nlayers), zomco(nlayers)
      real(kind=8) :: zrdnd(nlayers+1), zrdndc(nlayers+1)
      real(kind=8) :: zref(nlayers+1), zrefc(nlayers+1), zrefo(nlayers+1)
      real(kind=8) :: zrefd(nlayers+1), zrefdc(nlayers+1), zrefdo(nlayers+1)
      real(kind=8) :: zrup(nlayers+1), zrupd(nlayers+1)
      real(kind=8) :: zrupc(nlayers+1), zrupdc(nlayers+1)
      real(kind=8) :: zs1(nlayers+1)
      real(kind=8) :: ztauc(nlayers), ztauo(nlayers)
      real(kind=8) :: ztdn(nlayers+1), ztdnd(nlayers+1), ztdbt(nlayers+1)
      real(kind=8) :: ztoc(nlayers), ztor(nlayers)
      real(kind=8) :: ztra(nlayers+1), ztrac(nlayers+1), ztrao(nlayers+1)
      real(kind=8) :: ztrad(nlayers+1), ztradc(nlayers+1), ztrado(nlayers+1)
      real(kind=8) :: zdbtc(nlayers+1), ztdbtc(nlayers+1)
      real(kind=8) :: zincflx(ngptsw), zdbtc_nodel(nlayers+1) 
      real(kind=8) :: ztdbt_nodel(nlayers+1), ztdbtc_nodel(nlayers+1)

      real(kind=8) :: zdbtmc, zdbtmo, zf, zgw, zreflect
      real(kind=8) :: zwf, tauorig, repclc
!     real(kind=8) :: zincflux                                   ! inactive

! Arrays from rrtmg_sw_taumoln routines

!      real(kind=8) :: ztaug(nlayers,16), ztaur(nlayers,16)
!      real(kind=8) :: zsflxzen(16)
      real(kind=8) :: ztaug(nlayers,ngptsw), ztaur(nlayers,ngptsw)
      real(kind=8) :: zsflxzen(ngptsw)

! Arrays from rrtmg_sw_vrtqdr routine

      real(kind=8) :: zcd(nlayers+1,ngptsw), zcu(nlayers+1,ngptsw)
      real(kind=8) :: zfd(nlayers+1,ngptsw), zfu(nlayers+1,ngptsw)

! Inactive arrays
!     real(kind=8) :: zbbcd(nlayers+1), zbbcu(nlayers+1)
!     real(kind=8) :: zbbfd(nlayers+1), zbbfu(nlayers+1)
!     real(kind=8) :: zbbfddir(nlayers+1), zbbcddir(nlayers+1)

! ------------------------------------------------------------------

! Initializations

      ib1 = istart
      ib2 = iend
      klev = nlayers
      iw = 0
      repclc = 1.e-12_8
!      zincflux = 0.0_8

      do jk=1,klev+1
         pbbcd(jk)=0._8
         pbbcu(jk)=0._8
         pbbfd(jk)=0._8
         pbbfu(jk)=0._8
         pbbcddir(jk)=0._8
         pbbfddir(jk)=0._8
         puvcd(jk)=0._8
         puvfd(jk)=0._8
         puvcddir(jk)=0._8
         puvfddir(jk)=0._8
         pnicd(jk)=0._8
         pnifd(jk)=0._8
         pnicddir(jk)=0._8
         pnifddir(jk)=0._8
      enddo


! Calculate the optical depths for gaseous absorption and Rayleigh scattering

      call taumol_sw(klev, &
                     colh2o, colco2, colch4, colo2, colo3, colmol, &
                     laytrop, jp, jt, jt1, &
                     fac00, fac01, fac10, fac11, &
                     selffac, selffrac, indself, forfac, forfrac, indfor, &
                     zsflxzen, ztaug, ztaur)


! Top of shortwave spectral band loop, jb = 16 -> 29; ibm = 1 -> 14

      do jb = ib1, ib2
         ibm = jb-15
         igt = ngc(ibm)

! Reinitialize g-point counter for each band if output for each band is requested.
         if (iout.gt.0.and.ibm.ge.2) iw = ngs(ibm-1)

!        do jk=1,klev+1
!           zbbcd(jk)=0.0_8
!           zbbcu(jk)=0.0_8
!           zbbfd(jk)=0.0_8
!           zbbfu(jk)=0.0_8
!        enddo

! Top of g-point interval loop within each band (iw is cumulative counter) 
         do jg = 1,igt
            iw = iw+1

! Apply adjustments for correct Earth/Sun distance and zenith angle to incoming solar flux
            zincflx(iw) = adjflux(jb) * zsflxzen(iw) * prmu0
!             zincflux = zincflux + adjflux(jb) * zsflxzen(iw) * prmu0           ! inactive

! Compute layer reflectances and transmittances for direct and diffuse sources, 
! first clear then cloudy

! zrefc(jk)  direct albedo for clear
! zrefo(jk)  direct albedo for cloud
! zrefdc(jk) diffuse albedo for clear
! zrefdo(jk) diffuse albedo for cloud
! ztrac(jk)  direct transmittance for clear
! ztrao(jk)  direct transmittance for cloudy
! ztradc(jk) diffuse transmittance for clear
! ztrado(jk) diffuse transmittance for cloudy
!  
! zref(jk)   direct reflectance
! zrefd(jk)  diffuse reflectance
! ztra(jk)   direct transmittance
! ztrad(jk)  diffuse transmittance
!
! zdbtc(jk)  clear direct beam transmittance
! zdbto(jk)  cloudy direct beam transmittance
! zdbt(jk)   layer mean direct beam transmittance
! ztdbt(jk)  total direct beam transmittance at levels

! Clear-sky    
!   TOA direct beam    
            ztdbtc(1)=1.0_8
            ztdbtc_nodel(1)=1.0_8
!   Surface values
            zdbtc(klev+1) =0.0_8
            ztrac(klev+1) =0.0_8
            ztradc(klev+1)=0.0_8
            zrefc(klev+1) =palbp(ibm)
            zrefdc(klev+1)=palbd(ibm)
            zrupc(klev+1) =palbp(ibm)
            zrupdc(klev+1)=palbd(ibm)
           
! Cloudy-sky    
!   Surface values
            ztrao(klev+1) =0.0_8
            ztrado(klev+1)=0.0_8
            zrefo(klev+1) =palbp(ibm)
            zrefdo(klev+1)=palbd(ibm)
           
! Total sky    
!   TOA direct beam    
            ztdbt(1)=1.0_8
            ztdbt_nodel(1)=1.0_8
!   Surface values
            zdbt(klev+1) =0.0_8
            ztra(klev+1) =0.0_8
            ztrad(klev+1)=0.0_8
            zref(klev+1) =palbp(ibm)
            zrefd(klev+1)=palbd(ibm)
            zrup(klev+1) =palbp(ibm)
            zrupd(klev+1)=palbd(ibm)
    
    
! Top of layer loop
            do jk=1,klev

! Note: two-stream calculations proceed from top to bottom; 
!   RRTMG_SW quantities are given bottom to top and are reversed here

               ikl=klev+1-jk

! Set logical flag to do REFTRA calculation
!   Do REFTRA for all clear layers
               lrtchkclr(jk)=.true.

!   Do REFTRA only for cloudy layers in profile, since already done for clear layers
               lrtchkcld(jk)=.false.
               lrtchkcld(jk)=(pclfr(ikl) > repclc)

! Clear-sky optical parameters - this section inactive     
!   Original
!               ztauc(jk) = ztaur(ikl,iw) + ztaug(ikl,iw)
!               zomcc(jk) = ztaur(ikl,iw) / ztauc(jk)
!               zgcc(jk) = 0.0001_8
!   Total sky optical parameters        
!               ztauo(jk) = ztaur(ikl,iw) + ztaug(ikl,iw) + ptauc(ikl,ibm)
!               zomco(jk) = ptauc(ikl,ibm) * pomgc(ikl,ibm) + ztaur(ikl,iw)
!               zgco (jk) = (ptauc(ikl,ibm) * pomgc(ikl,ibm) * pasyc(ikl,ibm) + &
!                           ztaur(ikl,iw) * 0.0001_8) / zomco(jk)
!               zomco(jk) = zomco(jk) / ztauo(jk)

! Clear-sky optical parameters including aerosols
               ztauc(jk) = ztaur(ikl,iw) + ztaug(ikl,iw) + ptaua(ikl,ibm)
               zomcc(jk) = ztaur(ikl,iw) * 1.0_8 + ptaua(ikl,ibm) * pomga(ikl,ibm)
               zgcc(jk) = pasya(ikl,ibm) * pomga(ikl,ibm) * ptaua(ikl,ibm) / zomcc(jk)
               zomcc(jk) = zomcc(jk) / ztauc(jk)

! Pre-delta-scaling clear and cloudy direct beam transmittance (must use 'orig', unscaled cloud OD)       
!   \/\/\/ This block of code is only needed for unscaled direct beam calculation
               if (idelm .eq. 0) then
!     
                  zclear = 1.0_8 - pclfr(ikl)
                  zcloud = pclfr(ikl)

! Clear
!                   zdbtmc = exp(-ztauc(jk) / prmu0)
 
! Use exponential lookup table for transmittance, or expansion of exponential for low tau
                  ze1 = ztauc(jk) / prmu0
                  if (ze1 .le. od_lo) then
                     zdbtmc = 1._8 - ze1 + 0.5_8 * ze1 * ze1
                  else 
                     tblind = ze1 / (bpade + ze1)
                     itind = tblint * tblind + 0.5_8
                     zdbtmc = exp_tbl(itind)
                  endif

                  zdbtc_nodel(jk) = zdbtmc
                  ztdbtc_nodel(jk+1) = zdbtc_nodel(jk) * ztdbtc_nodel(jk)

! Clear + Cloud
                  tauorig = ztauc(jk) + ptaucorig(ikl,ibm)
!                   zdbtmo = exp(-tauorig / prmu0)

! Use exponential lookup table for transmittance, or expansion of exponential for low tau
                  ze1 = tauorig / prmu0
                  if (ze1 .le. od_lo) then
                     zdbtmo = 1._8 - ze1 + 0.5_8 * ze1 * ze1
                  else
                     tblind = ze1 / (bpade + ze1)
                     itind = tblint * tblind + 0.5_8
                     zdbtmo = exp_tbl(itind)
                  endif

                  zdbt_nodel(jk) = zclear * zdbtmc + zcloud * zdbtmo
                  ztdbt_nodel(jk+1) = zdbt_nodel(jk) * ztdbt_nodel(jk)

               endif
!   /\/\/\ Above code only needed for unscaled direct beam calculation


! Delta scaling - clear   
               zf = zgcc(jk) * zgcc(jk)
               zwf = zomcc(jk) * zf
               ztauc(jk) = (1.0_8 - zwf) * ztauc(jk)
               zomcc(jk) = (zomcc(jk) - zwf) / (1.0_8 - zwf)
               zgcc (jk) = (zgcc(jk) - zf) / (1.0_8 - zf)

! Total sky optical parameters (cloud properties already delta-scaled)
!   Use this code if cloud properties are derived in rrtmg_sw_cldprop       
               if (icpr .ge. 1) then
                  ztauo(jk) = ztauc(jk) + ptauc(ikl,ibm)
                  zomco(jk) = ztauc(jk) * zomcc(jk) + ptauc(ikl,ibm) * pomgc(ikl,ibm) 
                  zgco (jk) = (ptauc(ikl,ibm) * pomgc(ikl,ibm) * pasyc(ikl,ibm) + &
                              ztauc(jk) * zomcc(jk) * zgcc(jk)) / zomco(jk)
                  zomco(jk) = zomco(jk) / ztauo(jk)

! Total sky optical parameters (if cloud properties not delta scaled)
!   Use this code if cloud properties are not derived in rrtmg_sw_cldprop       
               elseif (icpr .eq. 0) then
                  ztauo(jk) = ztaur(ikl,iw) + ztaug(ikl,iw) + ptaua(ikl,ibm) + ptauc(ikl,ibm)
                  zomco(jk) = ptaua(ikl,ibm) * pomga(ikl,ibm) + ptauc(ikl,ibm) * pomgc(ikl,ibm) + &
                              ztaur(ikl,iw) * 1.0_8
                  zgco (jk) = (ptauc(ikl,ibm) * pomgc(ikl,ibm) * pasyc(ikl,ibm) + &
                              ptaua(ikl,ibm)*pomga(ikl,ibm)*pasya(ikl,ibm)) / zomco(jk)
                  zomco(jk) = zomco(jk) / ztauo(jk)

! Delta scaling - clouds 
!   Use only if subroutine rrtmg_sw_cldprop is not used to get cloud properties and to apply delta scaling
                  zf = zgco(jk) * zgco(jk)
                  zwf = zomco(jk) * zf
                  ztauo(jk) = (1._8 - zwf) * ztauo(jk)
                  zomco(jk) = (zomco(jk) - zwf) / (1.0_8 - zwf)
                  zgco (jk) = (zgco(jk) - zf) / (1.0_8 - zf)
               endif 

! End of layer loop
            enddo    


! Clear sky reflectivities
            call reftra_sw (klev, &
                            lrtchkclr, zgcc, prmu0, ztauc, zomcc, &
                            zrefc, zrefdc, ztrac, ztradc)

! Total sky reflectivities      
            call reftra_sw (klev, &
                            lrtchkcld, zgco, prmu0, ztauo, zomco, &
                            zrefo, zrefdo, ztrao, ztrado)


            do jk=1,klev

! Combine clear and cloudy contributions for total sky
               ikl = klev+1-jk 
               zclear = 1.0_8 - pclfr(ikl)
               zcloud = pclfr(ikl)

               zref(jk) = zclear*zrefc(jk) + zcloud*zrefo(jk)
               zrefd(jk)= zclear*zrefdc(jk) + zcloud*zrefdo(jk)
               ztra(jk) = zclear*ztrac(jk) + zcloud*ztrao(jk)
               ztrad(jk)= zclear*ztradc(jk) + zcloud*ztrado(jk)

! Direct beam transmittance        

! Clear
!                zdbtmc = exp(-ztauc(jk) / prmu0)

! Use exponential lookup table for transmittance, or expansion of 
! exponential for low tau
               ze1 = ztauc(jk) / prmu0
               if (ze1 .le. od_lo) then
                  zdbtmc = 1._8 - ze1 + 0.5_8 * ze1 * ze1
               else
                  tblind = ze1 / (bpade + ze1)
                  itind = tblint * tblind + 0.5_8
                  zdbtmc = exp_tbl(itind)
               endif

               zdbtc(jk) = zdbtmc
               ztdbtc(jk+1) = zdbtc(jk)*ztdbtc(jk)

! Clear + Cloud
!                zdbtmo = exp(-ztauo(jk) / prmu0)

! Use exponential lookup table for transmittance, or expansion of 
! exponential for low tau
               ze1 = ztauo(jk) / prmu0
               if (ze1 .le. od_lo) then
                  zdbtmo = 1._8 - ze1 + 0.5_8 * ze1 * ze1
               else
                  tblind = ze1 / (bpade + ze1)
                  itind = tblint * tblind + 0.5_8
                  zdbtmo = exp_tbl(itind)
               endif

               zdbt(jk) = zclear*zdbtmc + zcloud*zdbtmo
               ztdbt(jk+1) = zdbt(jk)*ztdbt(jk)
        
            enddo           
                 
! Vertical quadrature for clear-sky fluxes

            call vrtqdr_sw (klev, iw, &
                            zrefc, zrefdc, ztrac, ztradc, &
                            zdbtc, zrdndc, zrupc, zrupdc, ztdbtc, &
                            zcd, zcu)
      
! Vertical quadrature for cloudy fluxes

            call vrtqdr_sw (klev, iw, &
                            zref, zrefd, ztra, ztrad, &
                            zdbt, zrdnd, zrup, zrupd, ztdbt, &
                            zfd, zfu)

! Upwelling and downwelling fluxes at levels
!   Two-stream calculations go from top to bottom; 
!   layer indexing is reversed to go bottom to top for output arrays

            do jk=1,klev+1
               ikl=klev+2-jk

! Accumulate spectral fluxes over bands - inactive
!               zbbfu(ikl) = zbbfu(ikl) + zincflx(iw)*zfu(jk,iw)  
!               zbbfd(ikl) = zbbfd(ikl) + zincflx(iw)*zfd(jk,iw)
!               zbbcu(ikl) = zbbcu(ikl) + zincflx(iw)*zcu(jk,iw)
!               zbbcd(ikl) = zbbcd(ikl) + zincflx(iw)*zcd(jk,iw)
!               zbbfddir(ikl) = zbbfddir(ikl) + zincflx(iw)*ztdbt_nodel(jk)
!               zbbcddir(ikl) = zbbcddir(ikl) + zincflx(iw)*ztdbtc_nodel(jk)

! Accumulate spectral fluxes over whole spectrum  
               pbbfu(ikl) = pbbfu(ikl) + zincflx(iw)*zfu(jk,iw)
               pbbfd(ikl) = pbbfd(ikl) + zincflx(iw)*zfd(jk,iw)
               pbbcu(ikl) = pbbcu(ikl) + zincflx(iw)*zcu(jk,iw)
               pbbcd(ikl) = pbbcd(ikl) + zincflx(iw)*zcd(jk,iw)
               if (idelm .eq. 0) then 
                  pbbfddir(ikl) = pbbfddir(ikl) + zincflx(iw)*ztdbt_nodel(jk)
                  pbbcddir(ikl) = pbbcddir(ikl) + zincflx(iw)*ztdbtc_nodel(jk)
               elseif (idelm .eq. 1) then
                  pbbfddir(ikl) = pbbfddir(ikl) + zincflx(iw)*ztdbt(jk)
                  pbbcddir(ikl) = pbbcddir(ikl) + zincflx(iw)*ztdbtc(jk)
               endif

! Accumulate direct fluxes for UV/visible bands
               if (ibm >= 10 .and. ibm <= 13) then
                  puvcd(ikl) = puvcd(ikl) + zincflx(iw)*zcd(jk,iw)
                  puvfd(ikl) = puvfd(ikl) + zincflx(iw)*zfd(jk,iw)
                  if (idelm .eq. 0) then 
                     puvfddir(ikl) = puvfddir(ikl) + zincflx(iw)*ztdbt_nodel(jk)
                     puvcddir(ikl) = puvcddir(ikl) + zincflx(iw)*ztdbtc_nodel(jk)
                  elseif (idelm .eq. 1) then
                     puvfddir(ikl) = puvfddir(ikl) + zincflx(iw)*ztdbt(jk)
                     puvcddir(ikl) = puvcddir(ikl) + zincflx(iw)*ztdbtc(jk)
                  endif
! Accumulate direct fluxes for near-IR bands
               else if (ibm == 14 .or. ibm <= 9) then  
                  pnicd(ikl) = pnicd(ikl) + zincflx(iw)*zcd(jk,iw)
                  pnifd(ikl) = pnifd(ikl) + zincflx(iw)*zfd(jk,iw)
                  if (idelm .eq. 0) then 
                     pnifddir(ikl) = pnifddir(ikl) + zincflx(iw)*ztdbt_nodel(jk)
                     pnicddir(ikl) = pnicddir(ikl) + zincflx(iw)*ztdbtc_nodel(jk)
                  elseif (idelm .eq. 1) then
                     pnifddir(ikl) = pnifddir(ikl) + zincflx(iw)*ztdbt(jk)
                     pnicddir(ikl) = pnicddir(ikl) + zincflx(iw)*ztdbtc(jk)
                  endif
               endif

            enddo

! End loop on jg, g-point interval
         enddo             

! End loop on jb, spectral band
      enddo            
      
      end subroutine spcvrt_sw

      end module rrtmg_sw_spcvrt


