      module rrlw_kg10

      

      implicit none
      save

!-----------------------------------------------------------------
! rrtmg_lw ORIGINAL abs. coefficients for interval 10
! band 10:  1390-1480 cm-1 (low - h2o; high - h2o)
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
!-----------------------------------------------------------------
!
!  name     type     purpose
!  ----   : ----   : ---------------------------------------------
!fracrefao: real    
!fracrefbo: real    
! kao     : real     
! kbo     : real     
! selfrefo: real     
! forrefo : real     
!-----------------------------------------------------------------

      integer(kind=4), parameter :: no10 = 16

      real(kind=8) , dimension(no10) :: fracrefao
      real(kind=8) , dimension(no10) :: fracrefbo

      real(kind=8) :: kao(5,13,no10)
      real(kind=8) :: kbo(5,13:59,no10)
      real(kind=8) :: selfrefo(10,no10)
      real(kind=8) :: forrefo(4,no10)

!-----------------------------------------------------------------
! rrtmg_lw COMBINED abs. coefficients for interval 10
! band 10:  1390-1480 cm-1 (low - h2o; high - h2o)
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
!-----------------------------------------------------------------
!
!  name     type     purpose
!  ----   : ----   : ---------------------------------------------
!fracrefao: real    
!fracrefbo: real    
! kao     : real     
! kbo     : real     
! selfref : real     
! forref  : real     
!
! absa    : real
! absb    : real
!-----------------------------------------------------------------

      integer(kind=4), parameter :: ng10 = 6

      real(kind=8) , dimension(ng10) :: fracrefa
      real(kind=8) , dimension(ng10) :: fracrefb

      real(kind=8) :: ka(5,13,ng10)   , absa(65,ng10)
      real(kind=8) :: kb(5,13:59,ng10), absb(235,ng10)
      real(kind=8) :: selfref(10,ng10)
      real(kind=8) :: forref(4,ng10)

      equivalence (ka(1,1,1),absa(1,1)),(kb(1,13,1),absb(1,1))

      end module rrlw_kg10
