      module rrlw_cld


      implicit none
      save

!------------------------------------------------------------------
! rrtmg_lw cloud property coefficients

! Revised: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
!------------------------------------------------------------------

!  name     type     purpose
! -----  :  ----   : ----------------------------------------------
! abscld1:  real   : 
! absice0:  real   : 
! absice1:  real   : 
! absice2:  real   : 
! absice3:  real   : 
! absliq0:  real   : 
! absliq1:  real   : 
!------------------------------------------------------------------

      real(kind=8) :: abscld1
      real(kind=8) , dimension(2) :: absice0
      real(kind=8) , dimension(2,5) :: absice1
      real(kind=8) , dimension(43,16) :: absice2
      real(kind=8) , dimension(46,16) :: absice3
      real(kind=8) :: absliq0
      real(kind=8) , dimension(58,16) :: absliq1

      end module rrlw_cld

