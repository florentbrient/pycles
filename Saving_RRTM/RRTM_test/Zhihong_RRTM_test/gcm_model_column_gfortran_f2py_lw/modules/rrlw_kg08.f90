      module rrlw_kg08

      

      implicit none
      save

!-----------------------------------------------------------------
! rrtmg_lw ORIGINAL abs. coefficients for interval 8
! band 8:  1080-1180 cm-1 (low (i.e.>~300mb) - h2o; high - o3)
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
! kao_mco2: real     
! kbo_mco2: real     
! kao_mn2o: real     
! kbo_mn2o: real     
! kao_mo3 : real     
! selfrefo: real     
! forrefo : real     
! cfc12o  : real     
!cfc22adjo: real     
!-----------------------------------------------------------------

      integer(kind=4), parameter :: no8  = 16

      real(kind=8) , dimension(no8) :: fracrefao
      real(kind=8) , dimension(no8) :: fracrefbo
      real(kind=8) , dimension(no8) :: cfc12o
      real(kind=8) , dimension(no8) :: cfc22adjo

      real(kind=8) :: kao(5,13,no8)
      real(kind=8) :: kao_mco2(19,no8)
      real(kind=8) :: kao_mn2o(19,no8)
      real(kind=8) :: kao_mo3(19,no8)
      real(kind=8) :: kbo(5,13:59,no8)
      real(kind=8) :: kbo_mco2(19,no8)
      real(kind=8) :: kbo_mn2o(19,no8)
      real(kind=8) :: selfrefo(10,no8)
      real(kind=8) :: forrefo(4,no8)

!-----------------------------------------------------------------
! rrtmg_lw COMBINED abs. coefficients for interval 8
! band 8:  1080-1180 cm-1 (low (i.e.>~300mb) - h2o; high - o3)
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
! ka_mco2 : real     
! kb_mco2 : real     
! ka_mn2o : real     
! kb_mn2o : real     
! ka_mo3  : real     
! selfref : real     
! forref  : real     
! cfc12   : real     
! cfc22adj: real     
!
! absa    : real
! absb    : real
!-----------------------------------------------------------------

      integer(kind=4), parameter :: ng8  = 8

      real(kind=8) , dimension(ng8) :: fracrefa
      real(kind=8) , dimension(ng8) :: fracrefb
      real(kind=8) , dimension(ng8) :: cfc12
      real(kind=8) , dimension(ng8) :: cfc22adj

      real(kind=8) :: ka(5,13,ng8)    ,absa(65,ng8)
      real(kind=8) :: kb(5,13:59,ng8) ,absb(235,ng8)
      real(kind=8) :: ka_mco2(19,ng8)
      real(kind=8) :: ka_mn2o(19,ng8)
      real(kind=8) :: ka_mo3(19,ng8)
      real(kind=8) :: kb_mco2(19,ng8)
      real(kind=8) :: kb_mn2o(19,ng8)
      real(kind=8) :: selfref(10,ng8)
      real(kind=8) :: forref(4,ng8)

      equivalence (ka(1,1,1),absa(1,1)),(kb(1,13,1),absb(1,1))

      end module rrlw_kg08

