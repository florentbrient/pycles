      module rrsw_kg29

      
      use parrrsw, only : ng29

      implicit none
      save

!-----------------------------------------------------------------
! rrtmg_sw ORIGINAL abs. coefficients for interval 29
! band 29:  820-2600 cm-1 (low - h2o; high - co2)
!
! Initial version:  JJMorcrette, ECMWF, oct1999
! Revised: MJIacono, AER, jul2006
! Revised: MJIacono, AER, aug2008
!-----------------------------------------------------------------
!
!  name     type     purpose
!  ----   : ----   : ---------------------------------------------
! kao     : real     
! kbo     : real     
! selfrefo: real     
! forrefo : real     
!sfluxrefo: real     
! absh2oo : real     
! absco2o : real     
!-----------------------------------------------------------------

      integer(kind=4), parameter :: no29 = 16

      real(kind=8) :: kao(5,13,no29)
      real(kind=8) :: kbo(5,13:59,no29)
      real(kind=8) :: selfrefo(10,no29), forrefo(4,no29)
      real(kind=8) :: sfluxrefo(no29)
      real(kind=8) :: absh2oo(no29), absco2o(no29)

      real(kind=8) :: rayl

!-----------------------------------------------------------------
! rrtmg_sw COMBINED abs. coefficients for interval 29
! band 29:  820-2600 cm-1 (low - h2o; high - co2)
!
! Initial version:  JJMorcrette, ECMWF, oct1999
! Revised: MJIacono, AER, jul2006
! Revised: MJIacono, AER, aug2008
!-----------------------------------------------------------------
!
!  name     type     purpose
!  ----   : ----   : ---------------------------------------------
! ka      : real     
! kb      : real     
! selfref : real     
! forref  : real     
! sfluxref: real     
! absh2o  : real     
! absco2  : real     
!-----------------------------------------------------------------

      real(kind=8) :: ka(5,13,ng29), absa(65,ng29)
      real(kind=8) :: kb(5,13:59,ng29), absb(235,ng29)
      real(kind=8) :: selfref(10,ng29), forref(4,ng29)
      real(kind=8) :: sfluxref(ng29)
      real(kind=8) :: absh2o(ng29), absco2(ng29)

      equivalence (ka(1,1,1),absa(1,1)), (kb(1,13,1),absb(1,1))

      end module rrsw_kg29

