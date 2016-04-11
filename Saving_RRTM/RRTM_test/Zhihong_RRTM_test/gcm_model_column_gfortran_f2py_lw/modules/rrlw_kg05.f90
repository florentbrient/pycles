      module rrlw_kg05

      

      implicit none
      save

!-----------------------------------------------------------------
! rrtmg_lw ORIGINAL abs. coefficients for interval 5
! band 5:  700-820 cm-1 (low - h2o,co2; high - o3,co2)
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
! kao_mo3 : real     
! selfrefo: real     
! forrefo : real     
! ccl4o   : real
!-----------------------------------------------------------------

      integer(kind=4), parameter :: no5  = 16

      real(kind=8) :: fracrefao(no5,9) ,fracrefbo(no5,5)
      real(kind=8) :: kao(9,5,13,no5)
      real(kind=8) :: kbo(5,5,13:59,no5)
      real(kind=8) :: kao_mo3(9,19,no5)
      real(kind=8) :: selfrefo(10,no5)
      real(kind=8) :: forrefo(4,no5)
      real(kind=8) :: ccl4o(no5)

!-----------------------------------------------------------------
! rrtmg_lw COMBINED abs. coefficients for interval 5
! band 5:  700-820 cm-1 (low - h2o,co2; high - o3,co2)
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
! ka_mo3  : real     
! selfref : real     
! forref  : real     
! ccl4    : real
!
! absa    : real
! absb    : real
!-----------------------------------------------------------------

      integer(kind=4), parameter :: ng5  = 16

      real(kind=8) :: fracrefa(ng5,9) ,fracrefb(ng5,5)
      real(kind=8) :: ka(9,5,13,ng5)   ,absa(585,ng5)
      real(kind=8) :: kb(5,5,13:59,ng5),absb(1175,ng5)
      real(kind=8) :: ka_mo3(9,19,ng5)
      real(kind=8) :: selfref(10,ng5)
      real(kind=8) :: forref(4,ng5)
      real(kind=8) :: ccl4(ng5)
      
      equivalence (ka(1,1,1,1),absa(1,1)),(kb(1,1,13,1),absb(1,1))

      end module rrlw_kg05

