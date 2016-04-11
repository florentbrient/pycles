!     path:      $Source$
!     author:    $Author: miacono $
!     revision:  $Revision: 23308 $
!     created:   $Date: 2013-12-27 17:23:51 -0500 (Fri, 27 Dec 2013) $

      module rrtmg_sw_cldprop

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

      
      use parrrsw, only : nbndsw, jpband, jpb1, jpb2
      use rrsw_cld, only : extliq1, ssaliq1, asyliq1, &
                           extice2, ssaice2, asyice2, &
                           extice3, ssaice3, asyice3, fdlice3, &
                           abari, bbari, cbari, dbari, ebari, fbari
      use rrsw_wvn, only : wavenum1, wavenum2
      use rrsw_vsn, only : hvrcld, hnamcld

      implicit none

      contains

! ----------------------------------------------------------------------------
      subroutine cldprop_sw(nlayers, inflag, iceflag, liqflag, cldfrac, &
                            tauc, ssac, asmc, fsfc, ciwp, clwp, rei, rel, &
                            taucldorig, taucloud, ssacloud, asmcloud)
! ----------------------------------------------------------------------------

! Purpose: Compute the cloud optical properties for each cloudy layer.
! Note: Only inflag = 0 and inflag=2/liqflag=1/iceflag=1,2,3 are available;
! (Hu & Stamnes, Ebert and Curry, Key, and Fu) are implemented.

! ------- Input -------

      integer(kind=4), intent(in) :: nlayers         ! total number of layers
      integer(kind=4), intent(in) :: inflag          ! see definitions
      integer(kind=4), intent(in) :: iceflag         ! see definitions
      integer(kind=4), intent(in) :: liqflag         ! see definitions

      real(kind=8), intent(in) :: cldfrac(nlayers)         ! cloud fraction
                                                      !    Dimensions: (nlayers)
      real(kind=8), intent(in) :: ciwp(nlayers)            ! cloud ice water path
                                                      !    Dimensions: (nlayers)
      real(kind=8), intent(in) :: clwp(nlayers)            ! cloud liquid water path
                                                      !    Dimensions: (nlayers)
      real(kind=8), intent(in) :: rei(nlayers)             ! cloud ice particle effective size (microns)
                                                      !    Dimensions: (nlayers)
                                                      ! specific definition of rei depends on setting of iceflag:
                                                      ! iceflag = 0: (inactive)
                                                      !              
                                                      ! iceflag = 1: ice effective radius, r_ec, (Ebert and Curry, 1992),
                                                      !              r_ec range is limited to 13.0 to 130.0 microns
                                                      ! iceflag = 2: ice effective radius, r_k, (Key, Streamer Ref. Manual, 1996)
                                                      !              r_k range is limited to 5.0 to 131.0 microns
                                                      ! iceflag = 3: generalized effective size, dge, (Fu, 1996),
                                                      !              dge range is limited to 5.0 to 140.0 microns
                                                      !              [dge = 1.0315 * r_ec]
      real(kind=8), intent(in) :: rel(nlayers)             ! cloud liquid particle effective radius (microns)
                                                      !    Dimensions: (nlayers)
      real(kind=8), intent(in) :: tauc(nbndsw,nlayers)          ! cloud optical depth
                                                      !    Dimensions: (nbndsw,nlayers)
      real(kind=8), intent(in) :: ssac(nbndsw,nlayers)          ! single scattering albedo
                                                      !    Dimensions: (nbndsw,nlayers)
      real(kind=8), intent(in) :: asmc(nbndsw,nlayers)          ! asymmetry parameter
                                                      !    Dimensions: (nbndsw,nlayers)
      real(kind=8), intent(in) :: fsfc(nbndsw,nlayers)          ! forward scattering fraction
                                                      !    Dimensions: (nbndsw,nlayers)

! ------- Output -------

      real(kind=8), intent(out) :: taucloud(nlayers,jpband)     ! cloud optical depth (delta scaled)
                                                      !    Dimensions: (nlayers,jpband)
      real(kind=8), intent(out) :: taucldorig(nlayers,jpband)   ! cloud optical depth (non-delta scaled)
                                                      !    Dimensions: (nlayers,jpband)
      real(kind=8), intent(out) :: ssacloud(nlayers,jpband)     ! single scattering albedo (delta scaled)
                                                      !    Dimensions: (nlayers,jpband)
      real(kind=8), intent(out) :: asmcloud(nlayers,jpband)     ! asymmetry parameter (delta scaled)
                                                      !    Dimensions: (nlayers,jpband)

! ------- Local -------

!      integer(kind=4) :: ncbands
      integer(kind=4) :: ib, ib1, ib2, lay, istr, index, icx

      real(kind=8), parameter :: eps = 1.e-06     ! epsilon
      real(kind=8), parameter :: cldmin = 1.e-20  ! minimum value for cloud quantities
      real(kind=8) :: cwp                            ! total cloud water path
      real(kind=8) :: radliq                         ! cloud liquid droplet radius (microns)
      real(kind=8) :: radice                         ! cloud ice effective size (microns)
      real(kind=8) :: factor
      real(kind=8) :: fint
      real(kind=8) :: tauctot(nlayers)               ! band integrated cloud optical depth

      real(kind=8) :: taucldorig_a, ssacloud_a, taucloud_a, ffp, ffp1, ffpssa
      real(kind=8) :: tauiceorig, scatice, ssaice, tauice, tauliqorig, scatliq, ssaliq, tauliq

      real(kind=8) :: fdelta(jpb1:jpb2)
      real(kind=8) :: extcoice(jpb1:jpb2), gice(jpb1:jpb2)
      real(kind=8) :: ssacoice(jpb1:jpb2), forwice(jpb1:jpb2)
      real(kind=8) :: extcoliq(jpb1:jpb2), gliq(jpb1:jpb2)
      real(kind=8) :: ssacoliq(jpb1:jpb2), forwliq(jpb1:jpb2)

! Initialize

      hvrcld = '$Revision: 23308 $'

!      ncbands = 29
      ib1 = jpb1
      ib2 = jpb2
      tauctot(:) = 0._8

      do lay = 1, nlayers
         do ib = ib1 , ib2
            taucldorig(lay,ib) = tauc(ib-15,lay)
            taucloud(lay,ib) = 0.0_8
            ssacloud(lay,ib) = 1.0_8
            asmcloud(lay,ib) = 0.0_8
            tauctot(lay) = tauctot(lay) + tauc(ib-15,lay)
         enddo
      enddo

! Main layer loop
      do lay = 1, nlayers

         cwp = ciwp(lay) + clwp(lay)
         if (cldfrac(lay) .ge. cldmin .and. &
            (cwp .ge. cldmin .or. tauctot(lay) .ge. cldmin)) then

! (inflag=0): Cloud optical properties input directly
! Cloud optical properties already defined in tauc, ssac, asmc are unscaled;
! Apply delta-M scaling here
            if (inflag .eq. 0) then

               do ib = ib1 , ib2
                  taucldorig_a = tauc(ib-15,lay)
                  ffp = fsfc(ib-15,lay)
                  ffp1 = 1.0_8 - ffp
                  ffpssa = 1.0_8 - ffp * ssac(ib-15,lay)
                  ssacloud_a = ffp1 * ssac(ib-15,lay) / ffpssa
                  taucloud_a = ffpssa * taucldorig_a

                  taucldorig(lay,ib) = taucldorig_a
                  ssacloud(lay,ib) = ssacloud_a
                  taucloud(lay,ib) = taucloud_a
                  asmcloud(lay,ib) = (asmc(ib-15,lay) - ffp) / (ffp1)
               enddo

! (inflag=2): Separate treatement of ice clouds and water clouds.
            elseif (inflag .eq. 2) then       
               radice = rei(lay)

! Calculation of absorption coefficients due to ice clouds.
               if (ciwp(lay) .eq. 0.0_8) then
                  do ib = ib1 , ib2
                     extcoice(ib) = 0.0_8
                     ssacoice(ib) = 0.0_8
                     gice(ib)     = 0.0_8
                     forwice(ib)  = 0.0_8
                  enddo

! (iceflag = 1): 
! Note: This option uses Ebert and Curry approach for all particle sizes similar to
! CAM3 implementation, though this is somewhat ineffective for large ice particles
               elseif (iceflag .eq. 1) then
                  if (radice .lt. 13.0_8 .or. radice .gt. 130._8) stop &
                     'ICE RADIUS OUT OF BOUNDS'
                  do ib = ib1, ib2
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
                     extcoice(ib) = abari(icx) + bbari(icx)/radice
                     ssacoice(ib) = 1._8 - cbari(icx) - dbari(icx) * radice
                     gice(ib) = ebari(icx) + fbari(icx) * radice

! Check to ensure upper limit of gice is within physical limits for large particles
                     if (gice(ib) .ge. 1.0_8) gice(ib) = 1.0_8 - eps
                     forwice(ib) = gice(ib)*gice(ib)
! Check to ensure all calculated quantities are within physical limits.
                     if (extcoice(ib) .lt. 0.0_8) stop 'ICE EXTINCTION LESS THAN 0.0'
                     if (ssacoice(ib) .gt. 1.0_8) stop 'ICE SSA GRTR THAN 1.0'
                     if (ssacoice(ib) .lt. 0.0_8) stop 'ICE SSA LESS THAN 0.0'
                     if (gice(ib) .gt. 1.0_8) stop 'ICE ASYM GRTR THAN 1.0'
                     if (gice(ib) .lt. 0.0_8) stop 'ICE ASYM LESS THAN 0.0'
                  enddo

! For iceflag=2 option, ice particle effective radius is limited to 5.0 to 131.0 microns

               elseif (iceflag .eq. 2) then
                  if (radice .lt. 5.0_8 .or. radice .gt. 131.0_8) stop 'ICE RADIUS OUT OF BOUNDS'
                  factor = (radice - 2._8)/3._8
                  index = int(factor)
                  if (index .eq. 43) index = 42
                  fint = factor - real(index,kind=8)
                  do ib = ib1, ib2
                     extcoice(ib) = extice2(index,ib) + fint * &
                                   (extice2(index+1,ib) -  extice2(index,ib))
                     ssacoice(ib) = ssaice2(index,ib) + fint * &
                                   (ssaice2(index+1,ib) -  ssaice2(index,ib))
                     gice(ib) = asyice2(index,ib) + fint * &
                                   (asyice2(index+1,ib) -  asyice2(index,ib))
                     forwice(ib) = gice(ib)*gice(ib)
! Check to ensure all calculated quantities are within physical limits.
                     if (extcoice(ib) .lt. 0.0_8) stop 'ICE EXTINCTION LESS THAN 0.0'
                     if (ssacoice(ib) .gt. 1.0_8) stop 'ICE SSA GRTR THAN 1.0'
                     if (ssacoice(ib) .lt. 0.0_8) stop 'ICE SSA LESS THAN 0.0'
                     if (gice(ib) .gt. 1.0_8) stop 'ICE ASYM GRTR THAN 1.0'
                     if (gice(ib) .lt. 0.0_8) stop 'ICE ASYM LESS THAN 0.0'
                  enddo

! For iceflag=3 option, ice particle generalized effective size is limited to 5.0 to 140.0 microns

               elseif (iceflag .eq. 3) then
                  if (radice .lt. 5.0_8 .or. radice .gt. 140.0_8) stop 'ICE GENERALIZED EFFECTIVE SIZE OUT OF BOUNDS'
                  factor = (radice - 2._8)/3._8
                  index = int(factor)
                  if (index .eq. 46) index = 45
                  fint = factor - real(index,kind=8)
                  do ib = ib1 , ib2
                     extcoice(ib) = extice3(index,ib) + fint * &
                                   (extice3(index+1,ib) - extice3(index,ib))
                     ssacoice(ib) = ssaice3(index,ib) + fint * &
                                   (ssaice3(index+1,ib) - ssaice3(index,ib))
                     gice(ib) = asyice3(index,ib) + fint * &
                               (asyice3(index+1,ib) - asyice3(index,ib))
                     fdelta(ib) = fdlice3(index,ib) + fint * &
                                 (fdlice3(index+1,ib) - fdlice3(index,ib))
                     if (fdelta(ib) .lt. 0.0_8) stop 'FDELTA LESS THAN 0.0'
                     if (fdelta(ib) .gt. 1.0_8) stop 'FDELTA GT THAN 1.0'                     
                     forwice(ib) = fdelta(ib) + 0.5_8 / ssacoice(ib)
! See Fu 1996 p. 2067 
                     if (forwice(ib) .gt. gice(ib)) forwice(ib) = gice(ib)
! Check to ensure all calculated quantities are within physical limits.
                     if (extcoice(ib) .lt. 0.0_8) stop 'ICE EXTINCTION LESS THAN 0.0'
                     if (ssacoice(ib) .gt. 1.0_8) stop 'ICE SSA GRTR THAN 1.0'
                     if (ssacoice(ib) .lt. 0.0_8) stop 'ICE SSA LESS THAN 0.0'
                     if (gice(ib) .gt. 1.0_8) stop 'ICE ASYM GRTR THAN 1.0'
                     if (gice(ib) .lt. 0.0_8) stop 'ICE ASYM LESS THAN 0.0'
                  enddo

               endif
                  
! Calculation of absorption coefficients due to water clouds.
                if (clwp(lay) .eq. 0.0_8) then
                   do ib = ib1 , ib2
                      extcoliq(ib) = 0.0_8
                      ssacoliq(ib) = 0.0_8
                      gliq(ib) = 0.0_8
                      forwliq(ib) = 0.0_8
                   enddo

                elseif (liqflag .eq. 1) then
                   radliq = rel(lay)
                   if (radliq .lt. 2.5_8 .or. radliq .gt. 60._8) stop &
                      'LIQUID EFFECTIVE RADIUS OUT OF BOUNDS'
                   index = int(radliq - 1.5_8)
                   if (index .eq. 0) index = 1
                   if (index .eq. 58) index = 57
                   fint = radliq - 1.5_8 - real(index,kind=8)
                   do ib = ib1 , ib2
                      extcoliq(ib) = extliq1(index,ib) + fint * &
                                    (extliq1(index+1,ib) - extliq1(index,ib))
                      ssacoliq(ib) = ssaliq1(index,ib) + fint * &
                                    (ssaliq1(index+1,ib) - ssaliq1(index,ib))
                      if (fint .lt. 0._8 .and. ssacoliq(ib) .gt. 1._8) &
                                     ssacoliq(ib) = ssaliq1(index,ib)
                      gliq(ib) = asyliq1(index,ib) + fint * &
                                (asyliq1(index+1,ib) - asyliq1(index,ib))
                      forwliq(ib) = gliq(ib)*gliq(ib)
! Check to ensure all calculated quantities are within physical limits.
                      if (extcoliq(ib) .lt. 0.0_8) stop 'LIQUID EXTINCTION LESS THAN 0.0'
                      if (ssacoliq(ib) .gt. 1.0_8) stop 'LIQUID SSA GRTR THAN 1.0'
                      if (ssacoliq(ib) .lt. 0.0_8) stop 'LIQUID SSA LESS THAN 0.0'
                      if (gliq(ib) .gt. 1.0_8) stop 'LIQUID ASYM GRTR THAN 1.0'
                      if (gliq(ib) .lt. 0.0_8) stop 'LIQUID ASYM LESS THAN 0.0'
                   enddo
                endif

                do ib = ib1 , ib2
                   tauliqorig = clwp(lay) * extcoliq(ib)
                   tauiceorig = ciwp(lay) * extcoice(ib)
                   taucldorig(lay,ib) = tauliqorig + tauiceorig

                   ssaliq = ssacoliq(ib) * (1.0_8 - forwliq(ib)) / &
                           (1.0_8 - forwliq(ib) * ssacoliq(ib))
                   tauliq = (1.0_8 - forwliq(ib) * ssacoliq(ib)) * tauliqorig
                   ssaice = ssacoice(ib) * (1.0_8 - forwice(ib)) / &
                           (1.0_8 - forwice(ib) * ssacoice(ib))
                   tauice = (1.0_8 - forwice(ib) * ssacoice(ib)) * tauiceorig

                   scatliq = ssaliq * tauliq
                   scatice = ssaice * tauice

                   taucloud(lay,ib) = tauliq + tauice

! Ensure non-zero taucmc and scatice
                   if (taucloud(lay,ib).eq.0.0_8) taucloud(lay,ib) = cldmin
                   if (scatice.eq.0.0_8) scatice = cldmin

                   ssacloud(lay,ib) = (scatliq + scatice) / taucloud(lay,ib)

                   if (iceflag .eq. 3) then
! In accordance with the 1996 Fu paper, equation A.3, 
! the moments for ice were calculated depending on whether using spheres
! or hexagonal ice crystals.
                      istr = 1
                      asmcloud(lay,ib) = (1.0_8/(scatliq+scatice)) * &
                         (scatliq*(gliq(ib)**istr - forwliq(ib)) / &
                         (1.0_8 - forwliq(ib)) + scatice * ((gice(ib)-forwice(ib)) / &
                         (1.0_8 - forwice(ib)))**istr)
                   else 
! This code is the standard method for delta-m scaling. 
                      istr = 1
                      asmcloud(lay,ib) = (scatliq *  &
                         (gliq(ib)**istr - forwliq(ib)) / &
                         (1.0_8 - forwliq(ib)) + scatice * (gice(ib)**istr - forwice(ib)) / &
                         (1.0_8 - forwice(ib)))/(scatliq + scatice)
                   endif 

                enddo

            endif

         endif

! End layer loop
      enddo

      end subroutine cldprop_sw

      end module rrtmg_sw_cldprop


