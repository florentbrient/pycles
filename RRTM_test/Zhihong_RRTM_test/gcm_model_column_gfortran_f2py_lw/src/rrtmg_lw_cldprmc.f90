!     path:      $Source: /storm/rc1/cvsroot/rc/rrtmg_lw/src/rrtmg_lw_cldprmc.f90,v $
!     author:    $Author: miacono $
!     revision:  $Revision: 1.9 $
!     created:   $Date: 2011/04/08 20:25:00 $
!
      module rrtmg_lw_cldprmc

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

      use parrrtm, only : ngptlw, nbndlw
      use rrlw_cld, only: abscld1, absliq0, absliq1, &
                          absice0, absice1, absice2, absice3
      use rrlw_wvn, only: ngb
      use rrlw_vsn, only: hvrclc, hnamclc

      implicit none

      contains

! ------------------------------------------------------------------------------
      subroutine cldprmc(nlayers, inflag, iceflag, liqflag, cldfmc, &
                         ciwpmc, clwpmc, reicmc, relqmc, ncbands, taucmc)
! ------------------------------------------------------------------------------

! Purpose:  Compute the cloud optical depth(s) for each cloudy layer.

! ------- Input -------

      integer(kind=4), intent(in) :: nlayers         ! total number of layers
      integer(kind=4), intent(in) :: inflag          ! see definitions
      integer(kind=4), intent(in) :: iceflag         ! see definitions
      integer(kind=4), intent(in) :: liqflag         ! see definitions

      real(kind=8), intent(in) :: cldfmc(ngptlw,nlayers)        ! cloud fraction [mcica]
                                                      !    Dimensions: (ngptlw,nlayers)
      real(kind=8), intent(in) :: ciwpmc(ngptlw,nlayers)        ! cloud ice water path [mcica]
                                                      !    Dimensions: (ngptlw,nlayers)
      real(kind=8), intent(in) :: clwpmc(ngptlw,nlayers)        ! cloud liquid water path [mcica]
                                                      !    Dimensions: (ngptlw,nlayers)
      real(kind=8), intent(in) :: relqmc(nlayers)          ! liquid particle effective radius (microns)
                                                      !    Dimensions: (nlayers)
      real(kind=8), intent(in) :: reicmc(nlayers)          ! ice particle effective radius (microns)
                                                      !    Dimensions: (nlayers)
                                                      ! specific definition of reicmc depends on setting of iceflag:
                                                      ! iceflag = 0: ice effective radius, r_ec, (Ebert and Curry, 1992),
                                                      !              r_ec must be >= 10.0 microns
                                                      ! iceflag = 1: ice effective radius, r_ec, (Ebert and Curry, 1992),
                                                      !              r_ec range is limited to 13.0 to 130.0 microns
                                                      ! iceflag = 2: ice effective radius, r_k, (Key, Streamer Ref. Manual, 1996)
                                                      !              r_k range is limited to 5.0 to 131.0 microns
                                                      ! iceflag = 3: generalized effective size, dge, (Fu, 1996),
                                                      !              dge range is limited to 5.0 to 140.0 microns
                                                      !              [dge = 1.0315 * r_ec]

! ------- Output -------

      integer(kind=4), intent(out) :: ncbands        ! number of cloud spectral bands
      real(kind=8), intent(inout) :: taucmc(ngptlw,nlayers)     ! cloud optical depth [mcica]
                                                      !    Dimensions: (ngptlw,nlayers)

! ------- Local -------

      integer(kind=4) :: lay                         ! Layer index
      integer(kind=4) :: ib                          ! spectral band index
      integer(kind=4) :: ig                          ! g-point interval index
      integer(kind=4) :: index 
      integer(kind=4) :: icb(nbndlw)

      real(kind=8) :: abscoice(ngptlw)               ! ice absorption coefficients
      real(kind=8) :: abscoliq(ngptlw)               ! liquid absorption coefficients
      real(kind=8) :: cwp                            ! cloud water path
      real(kind=8) :: radice                         ! cloud ice effective size (microns)
      real(kind=8) :: factor                         ! 
      real(kind=8) :: fint                           ! 
      real(kind=8) :: radliq                         ! cloud liquid droplet radius (microns)
      real(kind=8), parameter :: eps = 1.e-6      ! epsilon
      real(kind=8), parameter :: cldmin = 1.e-20  ! minimum value for cloud quantities

! ------- Definitions -------

!     Explanation of the method for each value of INFLAG.  Values of
!     0 or 1 for INFLAG do not distingish being liquid and ice clouds.
!     INFLAG = 2 does distinguish between liquid and ice clouds, and
!     requires further user input to specify the method to be used to 
!     compute the aborption due to each.
!     INFLAG = 0:  For each cloudy layer, the cloud fraction and (gray)
!                  optical depth are input.  
!     INFLAG = 1:  For each cloudy layer, the cloud fraction and cloud
!                  water path (g/m2) are input.  The (gray) cloud optical 
!                  depth is computed as in CCM2.
!     INFLAG = 2:  For each cloudy layer, the cloud fraction, cloud 
!                  water path (g/m2), and cloud ice fraction are input.
!       ICEFLAG = 0:  The ice effective radius (microns) is input and the
!                     optical depths due to ice clouds are computed as in CCM3.
!       ICEFLAG = 1:  The ice effective radius (microns) is input and the
!                     optical depths due to ice clouds are computed as in 
!                     Ebert and Curry, JGR, 97, 3831-3836 (1992).  The 
!                     spectral regions in this work have been matched with
!                     the spectral bands in RRTM to as great an extent 
!                     as possible:  
!                     E&C 1      IB = 5      RRTM bands 9-16
!                     E&C 2      IB = 4      RRTM bands 6-8
!                     E&C 3      IB = 3      RRTM bands 3-5
!                     E&C 4      IB = 2      RRTM band 2
!                     E&C 5      IB = 1      RRTM band 1
!       ICEFLAG = 2:  The ice effective radius (microns) is input and the
!                     optical properties due to ice clouds are computed from
!                     the optical properties stored in the RT code,
!                     STREAMER v3.0 (Reference: Key. J., Streamer 
!                     User's Guide, Cooperative Institute for
!                     Meteorological Satellite Studies, 2001, 96 pp.).
!                     Valid range of values for re are between 5.0 and
!                     131.0 micron.
!       ICEFLAG = 3: The ice generalized effective size (dge) is input
!                    and the optical properties, are calculated as in
!                    Q. Fu, J. Climate, (1998). Q. Fu provided high resolution
!                    tables which were appropriately averaged for the
!                    bands in RRTM_LW.  Linear interpolation is used to
!                    get the coefficients from the stored tables.
!                    Valid range of values for dge are between 5.0 and
!                    140.0 micron.
!       LIQFLAG = 0:  The optical depths due to water clouds are computed as
!                     in CCM3.
!       LIQFLAG = 1:  The water droplet effective radius (microns) is input 
!                     and the optical depths due to water clouds are computed 
!                     as in Hu and Stamnes, J., Clim., 6, 728-742, (1993).
!                     The values for absorption coefficients appropriate for
!                     the spectral bands in RRTM have been obtained for a 
!                     range of effective radii by an averaging procedure 
!                     based on the work of J. Pinto (private communication).
!                     Linear interpolation is used to get the absorption 
!                     coefficients for the input effective radius.

      data icb /1,2,3,3,3,4,4,4,5, 5, 5, 5, 5, 5, 5, 5/

      hvrclc = '$Revision: 1.9 $'

      ncbands = 1

! This initialization is done in rrtmg_lw_subcol.F90.
!      do lay = 1, nlayers
!         do ig = 1, ngptlw
!            taucmc(ig,lay) = 0.0_8
!         enddo
!      enddo

! Main layer loop
      do lay = 1, nlayers

        do ig = 1, ngptlw
          cwp = ciwpmc(ig,lay) + clwpmc(ig,lay)
          if (cldfmc(ig,lay) .ge. cldmin .and. &
             (cwp .ge. cldmin .or. taucmc(ig,lay) .ge. cldmin)) then

! Ice clouds and water clouds combined.
            if (inflag .eq. 0) then
! Cloud optical depth already defined in taucmc, return to main program
               return

            elseif(inflag .eq. 1) then 
                stop 'INFLAG = 1 OPTION NOT AVAILABLE WITH MCICA'
!               cwp = ciwpmc(ig,lay) + clwpmc(ig,lay)
!               taucmc(ig,lay) = abscld1 * cwp

! Separate treatement of ice clouds and water clouds.
            elseif(inflag .eq. 2) then
               radice = reicmc(lay)

! Calculation of absorption coefficients due to ice clouds.
               if (ciwpmc(ig,lay) .eq. 0.0_8) then
                  abscoice(ig) = 0.0_8

               elseif (iceflag .eq. 0) then
                  if (radice .lt. 10.0_8) stop 'ICE RADIUS TOO SMALL'
                  abscoice(ig) = absice0(1) + absice0(2)/radice

               elseif (iceflag .eq. 1) then
                  if (radice .lt. 13.0_8 .or. radice .gt. 130._8) stop &
                      'ICE RADIUS OUT OF BOUNDS'
                  ncbands = 5
                  ib = icb(ngb(ig))
                  abscoice(ig) = absice1(1,ib) + absice1(2,ib)/radice

! For iceflag=2 option, ice particle effective radius is limited to 5.0 to 131.0 microns

               elseif (iceflag .eq. 2) then
                  if (radice .lt. 5.0_8 .or. radice .gt. 131.0_8) stop 'ICE RADIUS OUT OF BOUNDS'
                     ncbands = 16
                     factor = (radice - 2._8)/3._8
                     index = int(factor)
                     if (index .eq. 43) index = 42
                     fint = factor - real(index)
                     ib = ngb(ig)
                     abscoice(ig) = &
                         absice2(index,ib) + fint * &
                         (absice2(index+1,ib) - (absice2(index,ib))) 
               
! For iceflag=3 option, ice particle generalized effective size is limited to 5.0 to 140.0 microns

               elseif (iceflag .eq. 3) then
                  if (radice .lt. 5.0_8 .or. radice .gt. 140.0_8) stop 'ICE GENERALIZED EFFECTIVE SIZE OUT OF BOUNDS'
                     ncbands = 16
                     factor = (radice - 2._8)/3._8
                     index = int(factor)
                     if (index .eq. 46) index = 45
                     fint = factor - real(index)
                     ib = ngb(ig)
                     abscoice(ig) = &
                         absice3(index,ib) + fint * &
                         (absice3(index+1,ib) - (absice3(index,ib)))
   
               endif
                  
! Calculation of absorption coefficients due to water clouds.
               if (clwpmc(ig,lay) .eq. 0.0_8) then
                  abscoliq(ig) = 0.0_8

               elseif (liqflag .eq. 0) then
                   abscoliq(ig) = absliq0

               elseif (liqflag .eq. 1) then
                  radliq = relqmc(lay)
                  if (radliq .lt. 2.5_8 .or. radliq .gt. 60._8) stop &
                       'LIQUID EFFECTIVE RADIUS OUT OF BOUNDS'
                  index = int(radliq - 1.5_8)
                  if (index .eq. 0) index = 1
                  if (index .eq. 58) index = 57
                  fint = radliq - 1.5_8 - real(index)
                  ib = ngb(ig)
                  abscoliq(ig) = &
                        absliq1(index,ib) + fint * &
                        (absliq1(index+1,ib) - (absliq1(index,ib)))
               endif

               taucmc(ig,lay) = ciwpmc(ig,lay) * abscoice(ig) + &
                                clwpmc(ig,lay) * abscoliq(ig)

            endif
         endif
         enddo
      enddo

      end subroutine cldprmc

      end module rrtmg_lw_cldprmc
