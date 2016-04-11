      module rrsw_tbl

      

      implicit none
      save

!------------------------------------------------------------------
! rrtmg_sw lookup table arrays

! Initial version: MJIacono, AER, may2007
! Revised: MJIacono, AER, aug2007
! Revised: MJIacono, AER, aug2008
!------------------------------------------------------------------

!  name     type     purpose
! -----  :  ----   : ----------------------------------------------
! ntbl   :  integer: Lookup table dimension
! tblint :  real   : Lookup table conversion factor
! tau_tbl:  real   : Clear-sky optical depth 
! exp_tbl:  real   : Exponential lookup table for transmittance
! od_lo  :  real   : Value of tau below which expansion is used
!                  : in place of lookup table
! pade   :  real   : Pade approximation constant
! bpade  :  real   : Inverse of Pade constant
!------------------------------------------------------------------

      integer(kind=4), parameter :: ntbl = 10000

      real(kind=8), parameter :: tblint = 10000.0

      real(kind=8), parameter :: od_lo = 0.06

      real(kind=8) :: tau_tbl
      real(kind=8) , dimension(0:ntbl) :: exp_tbl

      real(kind=8), parameter :: pade = 0.278
      real(kind=8) :: bpade

      end module rrsw_tbl

