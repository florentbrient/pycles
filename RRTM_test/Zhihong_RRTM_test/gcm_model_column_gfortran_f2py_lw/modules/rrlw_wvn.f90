      module rrlw_wvn

      use parrrtm, only : nbndlw, mg, ngptlw, maxinpx

      implicit none
      save

!------------------------------------------------------------------
! rrtmg_lw spectral information

! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
!------------------------------------------------------------------

!  name     type     purpose
! -----  :  ----   : ----------------------------------------------
! ng     :  integer: Number of original g-intervals in each spectral band
! nspa   :  integer: For the lower atmosphere, the number of reference
!                    atmospheres that are stored for each spectral band
!                    per pressure level and temperature.  Each of these
!                    atmospheres has different relative amounts of the 
!                    key species for the band (i.e. different binary
!                    species parameters).
! nspb   :  integer: Same as nspa for the upper atmosphere
!wavenum1:  real   : Spectral band lower boundary in wavenumbers
!wavenum2:  real   : Spectral band upper boundary in wavenumbers
! delwave:  real   : Spectral band width in wavenumbers
! totplnk:  real   : Integrated Planck value for each band; (band 16
!                    includes total from 2600 cm-1 to infinity)
!                    Used for calculation across total spectrum
!totplk16:  real   : Integrated Planck value for band 16 (2600-3250 cm-1)
!                    Used for calculation in band 16 only if 
!                    individual band output requested
!totplnkderiv: real: Integrated Planck function derivative with respect
!                    to temperature for each band; (band 16
!                    includes total from 2600 cm-1 to infinity)
!                    Used for calculation across total spectrum
!totplk16deriv:real: Integrated Planck function derivative with respect
!                    to temperature for band 16 (2600-3250 cm-1)
!                    Used for calculation in band 16 only if 
!                    individual band output requested
!
! ngc    :  integer: The number of new g-intervals in each band
! ngs    :  integer: The cumulative sum of new g-intervals for each band
! ngm    :  integer: The index of each new g-interval relative to the
!                    original 16 g-intervals in each band
! ngn    :  integer: The number of original g-intervals that are 
!                    combined to make each new g-intervals in each band
! ngb    :  integer: The band index for each new g-interval
! wt     :  real   : RRTM weights for the original 16 g-intervals
! rwgt   :  real   : Weights for combining original 16 g-intervals 
!                    (256 total) into reduced set of g-intervals 
!                    (140 total)
! nxmol  :  integer: Number of cross-section molecules
! ixindx :  integer: Flag for active cross-sections in calculation
!------------------------------------------------------------------

      integer(kind=4) :: ng(nbndlw)
      integer(kind=4) :: nspa(nbndlw)
      integer(kind=4) :: nspb(nbndlw)

      real(kind=8) :: wavenum1(nbndlw)
      real(kind=8) :: wavenum2(nbndlw)
      real(kind=8) :: delwave(nbndlw)

      real(kind=8) :: totplnk(181,nbndlw)
      real(kind=8) :: totplk16(181)

      real(kind=8) :: totplnkderiv(181,nbndlw)
      real(kind=8) :: totplk16deriv(181)

      integer(kind=4) :: ngc(nbndlw)
      integer(kind=4) :: ngs(nbndlw)
      integer(kind=4) :: ngn(ngptlw)
      integer(kind=4) :: ngb(ngptlw)
      integer(kind=4) :: ngm(nbndlw*mg)

      real(kind=8) :: wt(mg)
      real(kind=8) :: rwgt(nbndlw*mg)

      integer(kind=4) :: nxmol
      integer(kind=4) :: ixindx(maxinpx)

      end module rrlw_wvn
