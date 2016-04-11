      module rrsw_kg24

      
      use parrrsw, only : ng24

      implicit none
      save

!-----------------------------------------------------------------
! rrtmg_sw ORIGINAL abs. coefficients for interval 24
! band 24: 12850-16000 cm-1 (low - h2o,o2; high - o2)
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
! abso3ao : real     
! abso3bo : real     
! raylao  : real     
! raylbo  : real     
!-----------------------------------------------------------------

      integer(kind=4), parameter :: no24 = 16

      real(kind=8) :: kao(9,5,13,no24)
      real(kind=8) :: kbo(5,13:59,no24)
      real(kind=8) :: selfrefo(10,no24), forrefo(3,no24)
      real(kind=8) :: sfluxrefo(no24,9)
      real(kind=8) :: abso3ao(no24), abso3bo(no24)
      real(kind=8) :: raylao(no24,9), raylbo(no24)

!-----------------------------------------------------------------
! rrtmg_sw COMBINED abs. coefficients for interval 24
! band 24: 12850-16000 cm-1 (low - h2o,o2; high - o2)
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
! absa    : real
! absb    : real
! selfref : real     
! forref  : real
! sfluxref: real     
! abso3a  : real     
! abso3b  : real     
! rayla   : real     
! raylb   : real     
!-----------------------------------------------------------------

      real(kind=8) :: ka(9,5,13,ng24), absa(585,ng24)
      real(kind=8) :: kb(5,13:59,ng24), absb(235,ng24)
      real(kind=8) :: selfref(10,ng24), forref(3,ng24)
      real(kind=8) :: sfluxref(ng24,9)
      real(kind=8) :: abso3a(ng24), abso3b(ng24)
      real(kind=8) :: rayla(ng24,9), raylb(ng24)

      equivalence (ka(1,1,1,1),absa(1,1)), (kb(1,13,1),absb(1,1))

      end module rrsw_kg24

