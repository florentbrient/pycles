      module rrlw_kg09

      

      implicit none
      save

!-----------------------------------------------------------------
! rrtmg_lw ORIGINAL abs. coefficients for interval 9
! band 9:  1180-1390 cm-1 (low - h2o,ch4; high - ch4)
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
! kao_mn2o: real     
! kbo_mn2o: real     
! selfrefo: real     
! forrefo : real     
!-----------------------------------------------------------------

      integer(kind=4), parameter :: no9  = 16

      real(kind=8) , dimension(no9) :: fracrefbo

      real(kind=8) :: fracrefao(no9,9)
      real(kind=8) :: kao(9,5,13,no9)
      real(kind=8) :: kbo(5,13:59,no9)
      real(kind=8) :: kao_mn2o(9,19,no9)
      real(kind=8) :: kbo_mn2o(19,no9)
      real(kind=8) :: selfrefo(10,no9)
      real(kind=8) :: forrefo(4,no9)

!-----------------------------------------------------------------
! rrtmg_lw COMBINED abs. coefficients for interval 9
! band 9:  1180-1390 cm-1 (low - h2o,ch4; high - ch4)
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
! ka_mn2o : real     
! kb_mn2o : real     
! selfref : real     
! forref  : real     
!
! absa    : real
! absb    : real
!-----------------------------------------------------------------

      integer(kind=4), parameter :: ng9  = 12

      real(kind=8) , dimension(ng9) :: fracrefb
      real(kind=8) :: fracrefa(ng9,9)
      real(kind=8) :: ka(9,5,13,ng9) ,absa(585,ng9)
      real(kind=8) :: kb(5,13:59,ng9) ,absb(235,ng9)
      real(kind=8) :: ka_mn2o(9,19,ng9)
      real(kind=8) :: kb_mn2o(19,ng9)
      real(kind=8) :: selfref(10,ng9)
      real(kind=8) :: forref(4,ng9)

      equivalence (ka(1,1,1,1),absa(1,1)),(kb(1,13,1),absb(1,1))

      end module rrlw_kg09
