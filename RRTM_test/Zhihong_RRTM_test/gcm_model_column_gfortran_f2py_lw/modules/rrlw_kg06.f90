      module rrlw_kg06

      

      implicit none
      save

!-----------------------------------------------------------------
! rrtmg_lw ORIGINAL abs. coefficients for interval 6
! band 6:  820-980 cm-1 (low - h2o; high - nothing)
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
!-----------------------------------------------------------------
!
!  name     type     purpose
!  ----   : ----   : ---------------------------------------------
!fracrefao: real    
! kao     : real     
! kao_mco2: real     
! selfrefo: real     
! forrefo : real     
!cfc11adjo: real
! cfc12o  : real
!-----------------------------------------------------------------

      integer(kind=4), parameter :: no6  = 16

      real(kind=8) , dimension(no6) :: fracrefao
      real(kind=8) :: kao(5,13,no6)
      real(kind=8) :: kao_mco2(19,no6)
      real(kind=8) :: selfrefo(10,no6)
      real(kind=8) :: forrefo(4,no6)

      real(kind=8) , dimension(no6) :: cfc11adjo
      real(kind=8) , dimension(no6) :: cfc12o

!-----------------------------------------------------------------
! rrtmg_lw COMBINED abs. coefficients for interval 6
! band 6:  820-980 cm-1 (low - h2o; high - nothing)
!
! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
!-----------------------------------------------------------------
!
!  name     type     purpose
!  ----   : ----   : ---------------------------------------------
!fracrefa : real    
! ka      : real     
! ka_mco2 : real     
! selfref : real     
! forref  : real     
!cfc11adj : real
! cfc12   : real
!
! absa    : real
!-----------------------------------------------------------------

      integer(kind=4), parameter :: ng6  = 8

      real(kind=8) , dimension(ng6) :: fracrefa
      real(kind=8) :: ka(5,13,ng6),absa(65,ng6)
      real(kind=8) :: ka_mco2(19,ng6)
      real(kind=8) :: selfref(10,ng6)
      real(kind=8) :: forref(4,ng6)

      real(kind=8) , dimension(ng6) :: cfc11adj
      real(kind=8) , dimension(ng6) :: cfc12

      equivalence (ka(1,1,1),absa(1,1))

      end module rrlw_kg06
