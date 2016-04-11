      module rrlw_kg13

      

      implicit none
      save

!-----------------------------------------------------------------
! rrtmg_lw ORIGINAL abs. coefficients for interval 13
! band 13:  2080-2250 cm-1 (low - h2o,n2o; high - nothing)
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
! kao_mco : real     
! kbo_mo3 : real     
! selfrefo: real     
! forrefo : real     
!-----------------------------------------------------------------

      integer(kind=4), parameter :: no13 = 16

      real(kind=8) , dimension(no13) :: fracrefbo

      real(kind=8) :: fracrefao(no13,9)
      real(kind=8) :: kao(9,5,13,no13)
      real(kind=8) :: kao_mco2(9,19,no13)
      real(kind=8) :: kao_mco(9,19,no13)
      real(kind=8) :: kbo_mo3(19,no13)
      real(kind=8) :: selfrefo(10,no13)
      real(kind=8) :: forrefo(4,no13)

!-----------------------------------------------------------------
! rrtmg_lw COMBINED abs. coefficients for interval 13
! band 13:  2080-2250 cm-1 (low - h2o,n2o; high - nothing)
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
! ka_mco  : real     
! kb_mo3  : real     
! selfref : real     
! forref  : real     
!
! absa    : real
!-----------------------------------------------------------------

      integer(kind=4), parameter :: ng13 = 4

      real(kind=8) , dimension(ng13) :: fracrefb

      real(kind=8) :: fracrefa(ng13,9)
      real(kind=8) :: ka(9,5,13,ng13) ,absa(585,ng13)
      real(kind=8) :: ka_mco2(9,19,ng13)
      real(kind=8) :: ka_mco(9,19,ng13)
      real(kind=8) :: kb_mo3(19,ng13)
      real(kind=8) :: selfref(10,ng13)
      real(kind=8) :: forref(4,ng13)

      equivalence (ka(1,1,1,1),absa(1,1))

      end module rrlw_kg13
