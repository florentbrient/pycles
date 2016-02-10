!     path:      $Source$
!     author:    $Author: miacono $
!     revision:  $Revision: 23308 $
!     created:   $Date: 2013-12-27 17:23:51 -0500 (Fri, 27 Dec 2013) $

      module rrtmg_sw_cldprmc

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

      
      use parrrsw, only : ngptsw, jpband, jpb1, jpb2
      use rrsw_cld, only : extliq1, ssaliq1, asyliq1, &
                           extice2, ssaice2, asyice2, &
                           extice3, ssaice3, asyice3, fdlice3, &
                           abari, bbari, cbari, dbari, ebari, fbari
      use rrsw_wvn, only : wavenum1, wavenum2, ngb
      use rrsw_vsn, only : hvrclc, hnamclc

      implicit none

      contains

! ----------------------------------------------------------------------------
      subroutine cldprmc_sw(nlayers, inflag, iceflag, liqflag, cldfmc, &
                            ciwpmc, clwpmc, reicmc, relqmc, &
                            taormc, taucmc, ssacmc, asmcmc, fsfcmc)
! ----------------------------------------------------------------------------

! Purpose: Compute the cloud optical properties for each cloudy layer
! and g-point interval for use by the McICA method.  
! Note: Only inflag = 0 and inflag=2/liqflag=1/iceflag=1,2,3 are available;
! (Hu & Stamnes, Ebert and Curry, Key, and Fu) are implemented. 

! ------- Input -------

      integer(kind=4), intent(in) :: nlayers         ! total number of layers
      integer(kind=4), intent(in) :: inflag          ! see definitions
      integer(kind=4), intent(in) :: iceflag         ! see definitions
      integer(kind=4), intent(in) :: liqflag         ! see definitions

      real(kind=8), intent(in) :: cldfmc(ngptsw,nlayers)        ! cloud fraction [mcica]
                                                      !    Dimensions: (ngptsw,nlayers)
      real(kind=8), intent(in) :: ciwpmc(ngptsw,nlayers)        ! cloud ice water path [mcica]
                                                      !    Dimensions: (ngptsw,nlayers)
      real(kind=8), intent(in) :: clwpmc(ngptsw,nlayers)        ! cloud liquid water path [mcica]
                                                      !    Dimensions: (ngptsw,nlayers)
      real(kind=8), intent(in) :: relqmc(nlayers)          ! cloud liquid particle effective radius (microns)
                                                      !    Dimensions: (nlayers)
      real(kind=8), intent(in) :: reicmc(nlayers)          ! cloud ice particle effective radius (microns)
                                                      !    Dimensions: (nlayers)
                                                      ! specific definition of reicmc depends on setting of iceflag:
                                                      ! iceflag = 0: (inactive)
                                                      !              
                                                      ! iceflag = 1: ice effective radius, r_ec, (Ebert and Curry, 1992),
                                                      !              r_ec range is limited to 13.0 to 130.0 microns
                                                      ! iceflag = 2: ice effective radius, r_k, (Key, Streamer Ref. Manual, 1996)
                                                      !              r_k range is limited to 5.0 to 131.0 microns
                                                      ! iceflag = 3: generalized effective size, dge, (Fu, 1996),
                                                      !              dge range is limited to 5.0 to 140.0 microns
                                                      !              [dge = 1.0315 * r_ec]
      real(kind=8), intent(in) :: fsfcmc(ngptsw,nlayers)        ! cloud forward scattering fraction
                                                      !    Dimensions: (ngptsw,nlayers)

! ------- Output -------

      real(kind=8), intent(inout) :: taucmc(ngptsw,nlayers)     ! cloud optical depth (delta scaled)
                                                      !    Dimensions: (ngptsw,nlayers)
      real(kind=8), intent(inout) :: ssacmc(ngptsw,nlayers)     ! single scattering albedo (delta scaled)
                                                      !    Dimensions: (ngptsw,nlayers)
      real(kind=8), intent(inout) :: asmcmc(ngptsw,nlayers)     ! asymmetry parameter (delta scaled)
                                                      !    Dimensions: (ngptsw,nlayers)
      real(kind=8), intent(out) :: taormc(ngptsw,nlayers)       ! cloud optical depth (non-delta scaled)
                                                      !    Dimensions: (ngptsw,nlayers)

! ------- Local -------

!      integer(kind=4) :: ncbands
      integer(kind=4) :: ib, lay, istr, index, icx, ig

      real(kind=8), parameter :: eps = 1.e-06     ! epsilon
      real(kind=8), parameter :: cldmin = 1.e-20  ! minimum value for cloud quantities
      real(kind=8) :: cwp                            ! total cloud water path
      real(kind=8) :: radliq                         ! cloud liquid droplet radius (microns)
      real(kind=8) :: radice                         ! cloud ice effective size (microns)
      real(kind=8) :: factor
      real(kind=8) :: fint

      real(kind=8) :: taucldorig_a, taucloud_a, ssacloud_a, ffp, ffp1, ffpssa
      real(kind=8) :: tauiceorig, scatice, ssaice, tauice, tauliqorig, scatliq, ssaliq, tauliq

      real(kind=8) :: fdelta(ngptsw)
      real(kind=8) :: extcoice(ngptsw), gice(ngptsw)
      real(kind=8) :: ssacoice(ngptsw), forwice(ngptsw)
      real(kind=8) :: extcoliq(ngptsw), gliq(ngptsw)
      real(kind=8) :: ssacoliq(ngptsw), forwliq(ngptsw)

! Initialize

      hvrclc = '$Revision: 23308 $'

! Some of these initializations are done elsewhere
      do lay = 1, nlayers
         do ig = 1, ngptsw
            taormc(ig,lay) = taucmc(ig,lay)
!            taucmc(ig,lay) = 0.0_8
!            ssacmc(ig,lay) = 1.0_8
!            asmcmc(ig,lay) = 0.0_8
         enddo
      enddo

! Main layer loop
      do lay = 1, nlayers

! Main g-point interval loop
         do ig = 1, ngptsw 
            cwp = ciwpmc(ig,lay) + clwpmc(ig,lay)
            if (cldfmc(ig,lay) .ge. cldmin .and. &
               (cwp .ge. cldmin .or. taucmc(ig,lay) .ge. cldmin)) then

! (inflag=0): Cloud optical properties input directly
               if (inflag .eq. 0) then
! Cloud optical properties already defined in taucmc, ssacmc, asmcmc are unscaled;
! Apply delta-M scaling here (using Henyey-Greenstein approximation)
                  taucldorig_a = taucmc(ig,lay)
                  ffp = fsfcmc(ig,lay)
                  ffp1 = 1.0_8 - ffp
                  ffpssa = 1.0_8 - ffp * ssacmc(ig,lay)
                  ssacloud_a = ffp1 * ssacmc(ig,lay) / ffpssa
                  taucloud_a = ffpssa * taucldorig_a

                  taormc(ig,lay) = taucldorig_a
                  ssacmc(ig,lay) = ssacloud_a
                  taucmc(ig,lay) = taucloud_a
                  asmcmc(ig,lay) = (asmcmc(ig,lay) - ffp) / (ffp1)

               elseif (inflag .eq. 1) then 
                  stop 'INFLAG = 1 OPTION NOT AVAILABLE WITH MCICA'

! (inflag=2): Separate treatement of ice clouds and water clouds.
               elseif (inflag .eq. 2) then       
                  radice = reicmc(lay)

! Calculation of absorption coefficients due to ice clouds.
                  if (ciwpmc(ig,lay) .eq. 0.0_8) then
                     extcoice(ig) = 0.0_8
                     ssacoice(ig) = 0.0_8
                     gice(ig)     = 0.0_8
                     forwice(ig)  = 0.0_8

! (iceflag = 1): 
! Note: This option uses Ebert and Curry approach for all particle sizes similar to
! CAM3 implementation, though this is somewhat unjustified for large ice particles
                  elseif (iceflag .eq. 1) then
                     if (radice .lt. 13.0_8 .or. radice .gt. 130._8) stop &
                        'ICE RADIUS OUT OF BOUNDS'
                     ib = ngb(ig)
                     if (wavenum2(ib) .gt. 1.43e04_8) then
                        icx = 1
                     elseif (wavenum2(ib) .gt. 7.7e03_8) then
                        icx = 2
                     elseif (wavenum2(ib) .gt. 5.3e03_8) then
                        icx = 3
                     elseif (wavenum2(ib) .gt. 4.0e03_8) then
                        icx = 4
                     elseif (wavenum2(ib) .ge. 2.5e03_8) then
                        icx = 5
                     endif
                     extcoice(ig) = (abari(icx) + bbari(icx)/radice)
                     ssacoice(ig) = 1._8 - cbari(icx) - dbari(icx) * radice
                     gice(ig) = ebari(icx) + fbari(icx) * radice
! Check to ensure upper limit of gice is within physical limits for large particles
                     if (gice(ig).ge.1._8) gice(ig) = 1._8 - eps
                     forwice(ig) = gice(ig)*gice(ig)
! Check to ensure all calculated quantities are within physical limits.
                     if (extcoice(ig) .lt. 0.0_8) stop 'ICE EXTINCTION LESS THAN 0.0'
                     if (ssacoice(ig) .gt. 1.0_8) stop 'ICE SSA GRTR THAN 1.0'
                     if (ssacoice(ig) .lt. 0.0_8) stop 'ICE SSA LESS THAN 0.0'
                     if (gice(ig) .gt. 1.0_8) stop 'ICE ASYM GRTR THAN 1.0'
                     if (gice(ig) .lt. 0.0_8) stop 'ICE ASYM LESS THAN 0.0'

! For iceflag=2 option, ice particle effective radius is limited to 5.0 to 131.0 microns

                  elseif (iceflag .eq. 2) then
                     if (radice .lt. 5.0_8 .or. radice .gt. 131.0_8) stop 'ICE RADIUS OUT OF BOUNDS'
                     factor = (radice - 2._8)/3._8
                     index = int(factor)
                     if (index .eq. 43) index = 42
                     fint = factor - real(index,kind=8)
                     ib = ngb(ig)
                     extcoice(ig) = extice2(index,ib) + fint * &
                                   (extice2(index+1,ib) -  extice2(index,ib))
                     ssacoice(ig) = ssaice2(index,ib) + fint * &
                                   (ssaice2(index+1,ib) -  ssaice2(index,ib))
                     gice(ig) = asyice2(index,ib) + fint * &
                                   (asyice2(index+1,ib) -  asyice2(index,ib))
                     forwice(ig) = gice(ig)*gice(ig)
! Check to ensure all calculated quantities are within physical limits.
                     if (extcoice(ig) .lt. 0.0_8) stop 'ICE EXTINCTION LESS THAN 0.0'
                     if (ssacoice(ig) .gt. 1.0_8) stop 'ICE SSA GRTR THAN 1.0'
                     if (ssacoice(ig) .lt. 0.0_8) stop 'ICE SSA LESS THAN 0.0'
                     if (gice(ig) .gt. 1.0_8) stop 'ICE ASYM GRTR THAN 1.0'
                     if (gice(ig) .lt. 0.0_8) stop 'ICE ASYM LESS THAN 0.0'

! For iceflag=3 option, ice particle generalized effective size is limited to 5.0 to 140.0 microns

                  elseif (iceflag .eq. 3) then
                     if (radice .lt. 5.0_8 .or. radice .gt. 140.0_8) stop 'ICE GENERALIZED EFFECTIVE SIZE OUT OF BOUNDS'
                     factor = (radice - 2._8)/3._8
                     index = int(factor)
                     if (index .eq. 46) index = 45
                     fint = factor - real(index,kind=8)
                     ib = ngb(ig)
                     extcoice(ig) = extice3(index,ib) + fint * &
                                   (extice3(index+1,ib) - extice3(index,ib))
                     ssacoice(ig) = ssaice3(index,ib) + fint * &
                                   (ssaice3(index+1,ib) - ssaice3(index,ib))
                     gice(ig) = asyice3(index,ib) + fint * &
                               (asyice3(index+1,ib) - asyice3(index,ib))
                     fdelta(ig) = fdlice3(index,ib) + fint * &
                                 (fdlice3(index+1,ib) - fdlice3(index,ib))
                     if (fdelta(ig) .lt. 0.0_8) stop 'FDELTA LESS THAN 0.0'
                     if (fdelta(ig) .gt. 1.0_8) stop 'FDELTA GT THAN 1.0'
                     forwice(ig) = fdelta(ig) + 0.5_8 / ssacoice(ig)
! See Fu 1996 p. 2067 
                     if (forwice(ig) .gt. gice(ig)) forwice(ig) = gice(ig)
! Check to ensure all calculated quantities are within physical limits.  
                     if (extcoice(ig) .lt. 0.0_8) stop 'ICE EXTINCTION LESS THAN 0.0'
                     if (ssacoice(ig) .gt. 1.0_8) stop 'ICE SSA GRTR THAN 1.0'
                     if (ssacoice(ig) .lt. 0.0_8) stop 'ICE SSA LESS THAN 0.0'
                     if (gice(ig) .gt. 1.0_8) stop 'ICE ASYM GRTR THAN 1.0'
                     if (gice(ig) .lt. 0.0_8) stop 'ICE ASYM LESS THAN 0.0'

                  endif

! Calculation of absorption coefficients due to water clouds.
                  if (clwpmc(ig,lay) .eq. 0.0_8) then
                     extcoliq(ig) = 0.0_8
                     ssacoliq(ig) = 0.0_8
                     gliq(ig) = 0.0_8
                     forwliq(ig) = 0.0_8

                  elseif (liqflag .eq. 1) then
                     radliq = relqmc(lay)
                     if (radliq .lt. 2.5_8 .or. radliq .gt. 60._8) stop &
                        'LIQUID EFFECTIVE RADIUS OUT OF BOUNDS'
                     index = int(radliq - 1.5_8)
                     if (index .eq. 0) index = 1
                     if (index .eq. 58) index = 57
                     fint = radliq - 1.5_8 - real(index,kind=8)
                     ib = ngb(ig)
                     extcoliq(ig) = extliq1(index,ib) + fint * &
                                   (extliq1(index+1,ib) - extliq1(index,ib))
                     ssacoliq(ig) = ssaliq1(index,ib) + fint * &
                                   (ssaliq1(index+1,ib) - ssaliq1(index,ib))
                     if (fint .lt. 0._8 .and. ssacoliq(ig) .gt. 1._8) &
                                    ssacoliq(ig) = ssaliq1(index,ib)
                     gliq(ig) = asyliq1(index,ib) + fint * &
                               (asyliq1(index+1,ib) - asyliq1(index,ib))
                     forwliq(ig) = gliq(ig)*gliq(ig)
! Check to ensure all calculated quantities are within physical limits.
                     if (extcoliq(ig) .lt. 0.0_8) stop 'LIQUID EXTINCTION LESS THAN 0.0'
                     if (ssacoliq(ig) .gt. 1.0_8) stop 'LIQUID SSA GRTR THAN 1.0'
                     if (ssacoliq(ig) .lt. 0.0_8) stop 'LIQUID SSA LESS THAN 0.0'
                     if (gliq(ig) .gt. 1.0_8) stop 'LIQUID ASYM GRTR THAN 1.0'
                     if (gliq(ig) .lt. 0.0_8) stop 'LIQUID ASYM LESS THAN 0.0'
                  endif
   
                  tauliqorig = clwpmc(ig,lay) * extcoliq(ig)
                  tauiceorig = ciwpmc(ig,lay) * extcoice(ig)
                  taormc(ig,lay) = tauliqorig + tauiceorig

                  ssaliq = ssacoliq(ig) * (1._8 - forwliq(ig)) / &
                          (1._8 - forwliq(ig) * ssacoliq(ig))
                  tauliq = (1._8 - forwliq(ig) * ssacoliq(ig)) * tauliqorig
                  ssaice = ssacoice(ig) * (1._8 - forwice(ig)) / &
                          (1._8 - forwice(ig) * ssacoice(ig))
                  tauice = (1._8 - forwice(ig) * ssacoice(ig)) * tauiceorig

                  scatliq = ssaliq * tauliq
                  scatice = ssaice * tauice
                  taucmc(ig,lay) = tauliq + tauice

! Ensure non-zero taucmc and scatice
                  if(taucmc(ig,lay).eq.0.) taucmc(ig,lay) = cldmin
                  if(scatice.eq.0.) scatice = cldmin

                  ssacmc(ig,lay) = (scatliq + scatice) / taucmc(ig,lay)

                  if (iceflag .eq. 3) then
! In accordance with the 1996 Fu paper, equation A.3, 
! the moments for ice were calculated depending on whether using spheres
! or hexagonal ice crystals.
! Set asymetry parameter to first moment (istr=1)
                     istr = 1
                     asmcmc(ig,lay) = (1.0_8/(scatliq+scatice))* &
                        (scatliq*(gliq(ig)**istr - forwliq(ig)) / &
                        (1.0_8 - forwliq(ig)) + scatice * ((gice(ig)-forwice(ig))/ &
                        (1.0_8 - forwice(ig)))**istr)

                  else 
! This code is the standard method for delta-m scaling. 
! Set asymetry parameter to first moment (istr=1)
                     istr = 1
                     asmcmc(ig,lay) = (scatliq *  &
                        (gliq(ig)**istr - forwliq(ig)) / &
                        (1.0_8 - forwliq(ig)) + scatice * (gice(ig)**istr - forwice(ig)) / &
                        (1.0_8 - forwice(ig)))/(scatliq + scatice)
                  endif 

               endif

            endif

! End g-point interval loop
         enddo

! End layer loop
      enddo

      end subroutine cldprmc_sw

      end module rrtmg_sw_cldprmc
