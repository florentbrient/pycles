!     path:      $Source$
!     author:    $Author: mike $
!     revision:  $Revision: 11661 $
!     created:   $Date: 2009-05-22 18:22:22 -0400 (Fri, 22 May 2009) $
!
      module rrtmg_sw_vrtqdr

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

      
      use parrrsw, only: ngptsw

      implicit none

      contains

! --------------------------------------------------------------------------
      subroutine vrtqdr_sw(klev, kw, &
                           pref, prefd, ptra, ptrad, &
                           pdbt, prdnd, prup, prupd, ptdbt, &
                           pfd, pfu)
! --------------------------------------------------------------------------
 
! Purpose: This routine performs the vertical quadrature integration
!
! Interface:  *vrtqdr_sw* is called from *spcvrt_sw* and *spcvmc_sw*
!
! Modifications.
! 
! Original: H. Barker
! Revision: Integrated with rrtmg_sw, J.-J. Morcrette, ECMWF, Oct 2002
! Revision: Reformatted for consistency with rrtmg_lw: MJIacono, AER, Jul 2006
!
!-----------------------------------------------------------------------

! ------- Declarations -------

! Input

      integer(kind=4), intent (in) :: klev                   ! number of model layers
      integer(kind=4), intent (in) :: kw                     ! g-point index

      real(kind=8), intent(in) :: pref(klev+1)                    ! direct beam reflectivity
                                                              !   Dimensions: (klev+1)
      real(kind=8), intent(in) :: prefd(klev+1)                   ! diffuse beam reflectivity
                                                              !   Dimensions: (klev+1)
      real(kind=8), intent(in) :: ptra(klev+1)                    ! direct beam transmissivity
                                                              !   Dimensions: (klev+1)
      real(kind=8), intent(in) :: ptrad(klev+1)                   ! diffuse beam transmissivity
                                                              !   Dimensions: (klev+1)

      real(kind=8), intent(in) :: pdbt(klev+1)
                                                              !   Dimensions: (klev+1)
      real(kind=8), intent(in) :: ptdbt(klev+1)
                                                              !   Dimensions: (klev+1)

      real(kind=8), intent(inout) :: prdnd(klev+1)
                                                              !   Dimensions: (klev+1)
      real(kind=8), intent(inout) :: prup(klev+1)
                                                              !   Dimensions: (klev+1)
      real(kind=8), intent(inout) :: prupd(klev+1)
                                                              !   Dimensions: (klev+1)

! Output
      real(kind=8), intent(out) :: pfd(klev+1,ngptsw)                  ! downwelling flux (W/m2)
                                                              !   Dimensions: (klev+1,ngptsw)
                                                              ! unadjusted for earth/sun distance or zenith angle
      real(kind=8), intent(out) :: pfu(klev+1,ngptsw)                  ! upwelling flux (W/m2)
                                                              !   Dimensions: (klev+1,ngptsw)
                                                              ! unadjusted for earth/sun distance or zenith angle

! Local

      integer(kind=4) :: ikp, ikx, jk

      real(kind=8) :: zreflect
      real(kind=8) :: ztdn(klev+1)  

! Definitions
!
! pref(jk)   direct reflectance
! prefd(jk)  diffuse reflectance
! ptra(jk)   direct transmittance
! ptrad(jk)  diffuse transmittance
!
! pdbt(jk)   layer mean direct beam transmittance
! ptdbt(jk)  total direct beam transmittance at levels
!
!-----------------------------------------------------------------------------
                   
! Link lowest layer with surface
             
      zreflect = 1._8 / (1._8 - prefd(klev+1) * prefd(klev))
      prup(klev) = pref(klev) + (ptrad(klev) * &
                 ((ptra(klev) - pdbt(klev)) * prefd(klev+1) + &
                   pdbt(klev) * pref(klev+1))) * zreflect
      prupd(klev) = prefd(klev) + ptrad(klev) * ptrad(klev) * &
                    prefd(klev+1) * zreflect

! Pass from bottom to top 

      do jk = 1,klev-1
         ikp = klev+1-jk                       
         ikx = ikp-1
         zreflect = 1._8 / (1._8 -prupd(ikp) * prefd(ikx))
         prup(ikx) = pref(ikx) + (ptrad(ikx) * &
                   ((ptra(ikx) - pdbt(ikx)) * prupd(ikp) + &
                     pdbt(ikx) * prup(ikp))) * zreflect
         prupd(ikx) = prefd(ikx) + ptrad(ikx) * ptrad(ikx) * &
                      prupd(ikp) * zreflect
      enddo
    
! Upper boundary conditions

      ztdn(1) = 1._8
      prdnd(1) = 0._8
      ztdn(2) = ptra(1)
      prdnd(2) = prefd(1)

! Pass from top to bottom

      do jk = 2,klev
         ikp = jk+1
         zreflect = 1._8 / (1._8 - prefd(jk) * prdnd(jk))
         ztdn(ikp) = ptdbt(jk) * ptra(jk) + &
                    (ptrad(jk) * ((ztdn(jk) - ptdbt(jk)) + &
                     ptdbt(jk) * pref(jk) * prdnd(jk))) * zreflect
         prdnd(ikp) = prefd(jk) + ptrad(jk) * ptrad(jk) * &
                      prdnd(jk) * zreflect
      enddo
    
! Up and down-welling fluxes at levels

      do jk = 1,klev+1
         zreflect = 1._8 / (1._8 - prdnd(jk) * prupd(jk))
         pfu(jk,kw) = (ptdbt(jk) * prup(jk) + &
                      (ztdn(jk) - ptdbt(jk)) * prupd(jk)) * zreflect
         pfd(jk,kw) = ptdbt(jk) + (ztdn(jk) - ptdbt(jk)+ &
                      ptdbt(jk) * prup(jk) * prdnd(jk)) * zreflect
      enddo

      end subroutine vrtqdr_sw

      end module rrtmg_sw_vrtqdr
