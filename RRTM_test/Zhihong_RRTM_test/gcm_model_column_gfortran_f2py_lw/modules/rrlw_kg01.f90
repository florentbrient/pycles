      module rrlw_kg01

      

      implicit none
      save

!-----------------------------------------------------------------
! rrtmg_lw ORIGINAL abs. coefficients for interval 1
! band 1:  10-250 cm-1 (low - h2o; high - h2o)
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
! kao_mn2 : real     
! kbo_mn2 : real     
! selfrefo: real     
! forrefo : real
!-----------------------------------------------------------------

      integer(kind=4), parameter :: no1  = 16

      real(kind=8) :: fracrefao(no1)  , fracrefbo(no1)
      real(kind=8) :: kao(5,13,no1)
      real(kind=8) :: kbo(5,13:59,no1)
      real(kind=8) :: kao_mn2(19,no1) , kbo_mn2(19,no1)
      real(kind=8) :: selfrefo(10,no1), forrefo(4,no1)

!-----------------------------------------------------------------
! rrtmg_lw COMBINED abs. coefficients for interval 1
! band 1:  10-250 cm-1 (low - h2o; high - h2o)
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
!-----------------------------------------------------------------
!
!  name     type     purpose
!  ----   : ----   : ---------------------------------------------
!fracrefa : real    
!fracrefb : real
! ka      : real     
! kb      : real     
! absa    : real
! absb    : real
! ka_mn2  : real     
! kb_mn2  : real     
! selfref : real     
! forref  : real
!-----------------------------------------------------------------

      integer(kind=4), parameter :: ng1  = 10

      real(kind=8) :: fracrefa(ng1)  , fracrefb(ng1)
      real(kind=8) :: ka(5,13,ng1)   , absa(65,ng1)
      real(kind=8) :: kb(5,13:59,ng1), absb(235,ng1)
      real(kind=8) :: ka_mn2(19,ng1) , kb_mn2(19,ng1)
      real(kind=8) :: selfref(10,ng1), forref(4,ng1)

      equivalence (ka(1,1,1),absa(1,1)), (kb(1,13,1),absb(1,1))

      end module rrlw_kg01




