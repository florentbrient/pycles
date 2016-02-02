module read_input

! -------- Modules --------

      implicit none
      public :: readprof, readcld, readaer

      contains
      
 !*************************************************************************
      real(kind=8) function earth_sun(idn)
!*************************************************************************
!
!  Purpose: Function to calculate the correction factor of Earth's orbit
!  for current day of the year

!  idn        : Day of the year
!  earth_sun  : square of the ratio of mean to actual Earth-Sun distance

! ------- Modules -------

      
      use rrsw_con, only : pi

      implicit none

      integer(kind=4), intent(in) :: idn

      real(kind=8) :: gamma

      gamma = 2._8*pi*(idn-1)/365._8

! Use Iqbal's equation 1.2.1

      earth_sun = 1.000110_8 + .034221_8 * cos(gamma) + .001289_8 * sin(gamma) + &
                   .000719_8 * cos(2._8*gamma) + .000077_8 * sin(2._8*gamma)

      end function earth_sun
           
!************************************************************************
      subroutine readprof(nlayers_out, iout_out, imca, icld_out, &
           iaer_out, isccos_out, idelm_out, pdp, &
           pavel_out, tavel_out, pz_out, tz_out, tbound_out, semiss_out, &
           zenith_out, adjflux_out, dyofyr_out, adjes_out, &
           coldry_out, wkl_out, inflag_out, iceflag_out,liqflag_out, &
           cldfrac_out, tauc, ssac, asmc, fsfc, ciwp, clwp, rei, rel, &
           tauaer_out, ssaaer_out, asmaer_out, filename, cldfile, aerfile)
!************************************************************************

! -------- Modules --------

      
      use rrsw_con, only: pi,planck,boltz,clight,avogad,alosmt,gascon, &
                          radcn1,radcn2
      use rrsw_wvn, only: wavenum1, wavenum2, delwave

! Note: COMMON blocks are left in this routine for array passing with
!       rrtatm.f for reading input profiles in single-column mode.
!       Scalars and arrays in subroutine call are renamed to avoid conflict
!       with COMMON blocks.
!                                                                         
! Purpose: Read in atmospheric profile.

      implicit integer(i-n), real(a-h,o-z)

! ------- Parameters -------
      parameter (mxlay = 203)
      parameter (nbndsw = 14)
      parameter (jpbands = 29)
      parameter (ib1 = 16, ib2 = 29)
      parameter (mg = 16)
      parameter (mxstr = 16)
      parameter (mcmu = 32)
      parameter (mxmol = 38)
      parameter (maxinpx = 35)
      parameter (maxxsec = 4)

      dimension altz(0:mxlay),ixtrans(14)
      dimension solvar(jpbands)

      common /control/  iaer, nstr, iout, istart, iend, icld, idelm, isccos
      common /constants/fluxfac,heatfac
      common /consts/   pic,planckc,boltzc,clightc,avogadc,alosmtc,gasconc, &
                        radcn1c,radcn2c
      common /swprop/   zenith, albedo, adjflux(jpbands)
      common /surface/  ireflect,semiss(jpbands)
      common /profile/  nlayers,pavel(mxlay),tavel(mxlay),pz(0:mxlay),tz(0:mxlay),tbound
      common /species/  coldry(mxlay),wkl(mxmol,mxlay),wbrodl(mxlay),colmol(mxlay),nmol
      common /ifil/     ird,ipr,ipu,idum(15)
      common /xsecctrl/ nxmol,ixindx(maxinpx)
      common /xsec/     wx(maxxsec,mxlay)
      common /pathx/    ixmax,nxmol0,ixindx0(maxinpx),wx0(maxinpx,mxlay)    
      common /xrrtatm/  ixsect

      common /cloudin/   inflag,clddat1(mxlay),clddat2(mxlay), &
                         iceflag,liqflag,clddat3(mxlay),clddat4(mxlay), &
                         clddatmom(0:16,mxlay)
      common /clouddat/  ncbands,cldfrac(mxlay), &
                         taucloud(mxlay,jpbands),ssacloud(mxlay,jpbands), &
                         xmom(0:16,mxlay,jpbands)
      common /aerdat/    ssaaer(mxlay,jpbands), phase(mcmu,mxlay,jpbands), &
                         tauaer(mxlay,jpbands)

      character*80 form1(0:1),form2(0:1),form3(0:1)
      character*1 ctest, cdollar, cdum

! Dimensions for transfer to rrtmg
      ! integer(kind=4), intent(in) :: ird_in              ! input file unit
      character(len=*), intent(in) :: filename           ! input file name
      character(len=*), intent(in), optional :: cldfile  ! input file name
      character(len=*), intent(in), optional :: aerfile  ! input file name
      integer(kind=4), intent(out) :: nlayers_out           ! total number of layers
      integer(kind=4), intent(out) :: imca                  ! McICA on/off flag (1 = use McICA)
      integer(kind=4), intent(out) :: icld_out              ! clear/cloud/overlap flag
      integer(kind=4), intent(out) :: iout_out              ! output option flag
      integer(kind=4), intent(out) :: iaer_out              ! aerosol flag
      integer(kind=4), intent(out) :: isccos_out            ! aerosol flag
      integer(kind=4), intent(out) :: idelm_out             ! aerosol flag

      real(kind=8), intent(out) :: pavel_out(mxlay)         ! layer pressures (mb) 
      real(kind=8), intent(out) :: tavel_out(mxlay)         ! layer temperatures (K)
      real(kind=8), intent(out) :: pz_out(0:mxlay)          ! level (interface) pressures (hPa, mb)
      real(kind=8), intent(out) :: tz_out(0:mxlay)          ! level (interface) temperatures (K)
      real(kind=8), intent(out) :: tbound_out               ! surface temperature (K)
      real(kind=8), intent(out) :: pdp(mxlay)               ! layer pressure thickness (hPa, mb)
      real(kind=8), intent(out) :: coldry_out(mxlay)        ! dry air molecular amount
      real(kind=8), intent(out) :: wkl_out(mxmol,mxlay)     ! molecular amounts (mol/cm-2)
      real(kind=8), intent(out) :: semiss_out(jpbands)       ! surface emissivity
      real(kind=8), intent(out) :: zenith_out               ! cos solar zenith angle
      real(kind=8), intent(out) :: adjflux_out(jpbands)      ! adjustment for current Earth/Sun distance

      integer(kind=4), intent(out) :: inflag_out              ! cloud property option flag
      integer(kind=4), intent(out) :: iceflag_out             ! ice cloud property flag
      integer(kind=4), intent(out) :: liqflag_out             ! liquid cloud property flag
      
      ! Output ADDED by ZTAN
      integer(kind=4), intent(out) :: dyofyr_out
      real(kind=8), intent(out) :: adjes_out
      
      real(kind=8), intent(out) :: cldfrac_out(mxlay)         ! cloud fraction
      real(kind=8), intent(out) :: tauc(nbndsw,mxlay)         ! in-cloud optical depth (non-delta scaled)
      real(kind=8), intent(out) :: ssac(nbndsw,mxlay)         ! in-cloud single scattering albedo (non-delta scaled)
      real(kind=8), intent(out) :: asmc(nbndsw,mxlay)         ! in-cloud asymmetry parameter (non-delta scaled)
      real(kind=8), intent(out) :: fsfc(nbndsw,mxlay)         ! in-cloud forward scattering fraction (non-delta scaled)
      real(kind=8), intent(out) :: ciwp(mxlay)                ! in-cloud ice water path
      real(kind=8), intent(out) :: clwp(mxlay)                ! in-cloud liquid water path
      real(kind=8), intent(out) :: rei(mxlay)                 ! cloud ice particle size
      real(kind=8), intent(out) :: rel(mxlay)                 ! cloud liquid particle size
      real(kind=8), intent(out) :: tauaer_out(mxlay,jpbands)   ! aerosol optical depth
      real(kind=8), intent(out) :: ssaaer_out(mxlay,jpbands)   ! aerosol single scattering albedo
      real(kind=8), intent(out) :: asmaer_out(mxlay,jpbands)   ! aerosol asymmetry parameter
                                                                 !   first momemnt of input phase function
!
! Local
      integer(kind=4) :: juldat                               ! day of year
      real(kind=8) :: fice(mxlay)                             ! cloud ice fraction

! Initializations

      data cdollar /'$'/
      data ixtrans /0,0,0,1,2,3,0,0,0,0,0,4,0,0/

      pi = 2._8 * asin(1._8)   ! Added by ZTAN 
      
      form1(0) = '(3f10.4,a3,i2,1x,2(f7.2,f8.3,f7.2))'
      form2(0) = '(3f10.4,a3,i2,23x,(f7.2,f8.3,f7.2))'
      form3(0) = '(8e10.3)'
      form1(1) = '(g15.7,g10.4,g10.4,a3,i2,1x,2(g7.2,g8.3,g7.2))'
      form2(1) = '(g15.7,g10.4,g10.4,a3,i2,23x,(g7.2,g8.3,g7.2))'
      form3(1) = '(8g15.7)'

! Pass constants to common block names for rrtatm
      pic = pi
      planckc = planck
      boltzc = boltz
      clightc = clight
      avogadc = avogad
      alosmtc = alosmt
      gasconc = gascon
      radcn1c = radcn1
      radcn2c = radcn2
 
      ixmax = maxinpx

! Open the input set of atmospheres
      open (ird,file=filename,form='formatted')

      do ilay = 1,mxlay
         do isp = 1,mxmol
            wkl(isp,ilay) = 0.0_8
         enddo
         do isp = 1,maxxsec
            wx(isp,ilay) = 0.0_8
         enddo
      enddo

! Top of read input loop
 1000 continue
      read (ird,9009,end=8800) ctest
      if (ctest .ne. cdollar) goto 1000

      read (ird,9011) iaer, iatm, iscat, istrm, iout, imca, icld, idelm, icos

      if (idelm.gt.1 .or. idelm.lt.0 .or. icos.gt.0 .or. icos.lt.0) then
         print *,'INVALID MEASUREMENT COMPARISON FLAG'
         stop
      endif
      isccos = icos

! No cross-sections implemented in shortwave.
      ixsect = 0

! Only 2-stream scattering implemented in rrtmg
      if (iscat .ne. 1) then
         print *,'INVALID SCATTERING OPTION CHOSEN'
         stop
      endif

      if (istrm .eq. 0) then 
         nstr = 2
      else 
         print *, 'INVALID VALUE FOR ISTRM'
         stop
      endif

      read (ird,9020) juldat, sza, isolvar, solvar(ib1:ib2)

      zenith = cos(sza * pi / 180._8)
      if (juldat .eq. 0) then
         adjflux_jd = 1._8
      else
         adjflux_jd = earth_sun (juldat)
      endif

! If clouds are present, read in appropriate input file, IN_CLD_RRTM.
      if (icld .ge. 1) call readcld(cldfile)

! If aerosols are present, read in appropriate input from file, IN_AER_RRTM. 
      if (iaer .eq. 10) call readaer(aerfile)


      if (isolvar .eq. 0) then
         do ib = ib1,ib2
            adjflux(ib) = adjflux_jd
         enddo
      elseif (isolvar .eq. 1) then
         do ib=ib1,ib2
            adjflux(ib) = adjflux_jd * solvar(ib1)
         enddo
      elseif (isolvar .eq. 2) then
         do ib=ib1,ib2
            adjflux(ib) = adjflux_jd * solvar(ib)
         enddo
      else
         print *, 'ISOLVAR = ', isolvar, ' NOT A VALID INPUT VALUE'
         stop
      endif
      
      ! Output added by ZTAN:
      dyofyr_out = juldat
      adjes_out  = adjflux_jd

      read (ird,9012) iemis, ireflect, semiss(ib1:ib2)
      if (iemis .eq. 0) then
         do ib = ib1, ib2
            semiss(ib) = 1._8
         enddo
      elseif (iemis .eq. 1) then
         do ib = ib1, ib2
            semiss(ib) = semiss(ib1)
         enddo
      elseif (iemis .eq. 2) then
!          print *, 'THESE ARE THE INPUT EMISSIVITY VALUES'
!          print *, semiss(ib1:ib2)
      else
          print *, 'IEMIS = ', iemis, ' NOT A VALID INPUT VALUE'
          stop
      endif
     
      if (iatm .eq. 0) then
         read (ird,9013) iform,nlayers,nmol
         if (nmol.eq.0) nmol = 7                                    
         read (ird,form1(iform)) pavel(1),tavel(1),secntk,cinp, &
              ipthak,altz(0),pz(0),tz(0),altz(1),pz(1),tz(1)
         read (ird,form3(iform)) (wkl(m,1),m=1,7), wbrodl(1)
         if(nmol .gt. 7) read (ird,form3(iform)) (wkl(m,1),m=8,nmol)

         do l = 2, nlayers
            read (ird,form2(iform)) pavel(l),tavel(l),secntk,cinp, &
                 ipthrk,altz(l),pz(l),tz(l)
            read (ird,form3(iform)) (wkl(m,l),m=1,7), wbrodl(l)
            if(nmol .gt. 7) read (ird,form3(iform)) (wkl(m,l),m=8,nmol)
         enddo
           
         if (ixsect .eq. 1) then                                 
            read (ird,9300) nxmol0
            nxmol = nxmol0
            call xsident(ird)
            read (ird,9301) iformx
     
            do l = 1, nlayers       
               read (ird,9010) cdum
               read (ird, form3(iformx)) (wx0(m,l),m=1,7),wbrodx    
               if (nxmol0 .gt. 7) read (ird,form3(iformx)) &
                  (wx0(m,l),m=8,nxmol0)
            enddo
         endif
      else
         ipu = 7
         ipr = 66
         open(unit=ipr,file='tape6',status='unknown')
         call rrtatm
         if (ixsect .eq. 1) then
            do mx = 1, nxmol0
               ixindx(mx) = ixtrans(ixindx0(mx))
            enddo
         endif
      endif

! Test for mixing ratio input.
      imix = 1
      do m = 1, nmol
         if (wkl(m,1) .gt. 1.0_8) then
            imix = 0
            goto 3600
         endif
      enddo
 3600 continue

      if (ixsect .eq. 1) then
         imixx = 0
         if (wx0(1,1) .le. 1.0_8) imixx = 1
      endif
      do l = 1, nlayers
         summol = 0.0_8
         do imol = 2, nmol
            summol = summol + wkl(imol,l)
         enddo
         if (imix .eq. 1) then
            coldry(l) = wbrodl(l) / (1._8 - summol)
            do imol = 1, nmol
               wkl(imol,l) = coldry(l) * wkl(imol,l)
            enddo
         else
            coldry(l) = wbrodl(l) + summol
         endif
         if (ixsect .eq. 1) then
            do ix = 1, nxmol0
               if (ixindx(ix) .ne. 0) then
                  if (imixx .eq. 1) then
                     wx(ixindx(ix),l) = coldry(l) * wx0(ix,l) * 1.e-20_8
                  else
                     wx(ixindx(ix),l) = wx0(ix,l) * 1.e-20_8
                  endif
               endif
            enddo 
         endif
      enddo

! Pass output arrays to new variables for transfer to rrtmg thorugh subroutine call.
      nlayers_out = nlayers
      iout_out = iout
      icld_out = icld
      iaer_out = iaer
      isccos_out = isccos
      idelm_out = idelm
      inflag_out = inflag
      iceflag_out = iceflag
      liqflag_out = liqflag

      pz_out(0) = pz(0)
      tz_out(0) = tz(0)
      tbound_out = tz(0)
      do l = 1, mxlay
         pavel_out(l) = pavel(l)
         tavel_out(l) = tavel(l)
         pz_out(l) = pz(l)
         tz_out(l) = tz(l)
         pdp(l) = (pz(l-1) - pz(l))
         coldry_out(l) = coldry(l)
         cldfrac_out(l) = cldfrac(l)
         do imol = 1,nmol
            wkl_out(imol,l) = wkl(imol,l)
         enddo
      enddo
      
      do l = 1, mxlay
         if (inflag.eq.0) then
            do n = 1, nbndsw
               tauc(n,l) = clddat1(l)
               ssac(n,l) = clddat2(l)
               asmc(n,l) = clddatmom(1,l)
               fsfc(n,l) = asmc(n,l)**2
            enddo
            ciwp(l) = 0._8
            clwp(l) = 0._8
            fice(l) = 0._8
            rei(l) = 0._8
            rel(l) = 0._8
         else
            do n = 1, nbndsw
               tauc(n,l) = 0._8
               ssac(n,l) = 1._8
               asmc(n,l) = 0._8
               fsfc(n,l) = 0._8
            enddo
            cwp = clddat1(l)
            fice(l) = clddat2(l)
            ciwp(l) = cwp * fice(l)
            clwp(l) = cwp * (1._8 - fice(l))
            rei(l) = clddat3(l)
            rel(l) = clddat4(l)
         endif 
      enddo

      do l = 1, mxlay
         do nb = ib1,ib2
            tauaer_out(l,nb) = tauaer(l,nb)
            ssaaer_out(l,nb) = ssaaer(l,nb)
            asmaer_out(l,nb) = phase(1,l,nb)
         enddo
      enddo
      zenith_out = zenith
      ! write(*,*) sza, pi, zenith, zenith_out
      do nb = ib1,ib2
        semiss_out(nb) = semiss(nb)
        adjflux_out(nb) = adjflux(nb)
      enddo
      
      goto 9000

 8800 continue
      stop ' INVALID INPUT_RRTM '

 9000 continue

 9009 format (a1,1x,i2,i2,i2)
 9010 format (a1)
 9011 format (18x,i2,29x,i1,32x,i1,1x,i1,2x,i3,3x,i1,i1,3x,i1,i1)
 9012 format (11x,i1,2x,i1,14f5.3)
 9013 format (1x,i1,i3,i5)                                     
 9020 format (12x, i3, 3x, f7.4, 4x, i1, 14f7.5)
 9300 format (i5)
 9301 format (1x,i1)

      end subroutine readprof

!***************************************************************************
      subroutine readcld(cldfile)
!***************************************************************************

      

! Purpose:  To read in IN_CLD_RRTM_SW, the file that contains input 
!           cloud properties.

      implicit integer(i-n), real(a-h,o-z)

      character(len=*), intent(in) :: cldfile           ! input file name
      
! ------- Parameters ------- 
      parameter (mxlay = 203, jpbands = 29)
      parameter (ib1 = 16, ib2 = 29)
      parameter (mg = 16)
      parameter (mxstr = 16)

      common /control/   iaer, nstr, iout, istart, iend, icld, idelm, isccos
      common /profile/   nlayers,pavel(mxlay),tavel(mxlay),pz(0:mxlay),tz(0:mxlay),tbound
      common /cloudin/   inflag,clddat1(mxlay),clddat2(mxlay), &
                         iceflag,liqflag,clddat3(mxlay),clddat4(mxlay), &
                         clddatmom(0:16,mxlay)
      common /clouddat/  ncbands,cldfrac(mxlay), &
                         taucloud(mxlay,jpbands),ssacloud(mxlay,jpbands), &
                         xmom(0:16,mxlay,jpbands)

      character*1 ctest, cpercent

      data cpercent /'%'/
      irdcld = 11

      open(irdcld,file=cldfile,form='formatted')

! Read in cloud input option.  
      read(irdcld,9050) inflag, iceflag, liqflag

      do lay = 1, nlayers
         cldfrac(lay) = 0._8
      enddo

      if (inflag .eq. 0) then
 950     continue
!  For INFLAG = 0 or 1, for each cloudy layer only LAY, FRAC, and
!  DAT1 are pertinent.  If CTEST = '%', then there are no more 
!  cloudy layers to process.
         read (irdcld,9099,end=8950) ctest,lay,frac,dat1,dat2,clddatmom(0:nstr,lay)
         if (ctest .eq. cpercent) goto 8950
         cldfrac(lay) = frac
         clddat1(lay) = dat1
         clddat2(lay) = dat2
         goto 950
 8950    continue

      else
 1000    continue
! For INFLAG = 0 or 1, for each cloudy layer only LAY, FRAC, and
! DAT1 are pertinent.  If CTEST = '%', then there are no more 
! cloudy layers to process.
         read (irdcld,9100,end=9000) ctest,lay,frac,dat1,dat2,dat3,dat4
         if (ctest .eq. cpercent) goto 9000
         cldfrac(lay) = frac
         clddat1(lay) = dat1
         clddat2(lay) = dat2
         clddat3(lay) = dat3
         clddat4(lay) = dat4
         goto 1000
 9000    continue
      endif

      close(irdcld)

 9050 format (3x,i2,4x,i1,4x,i1)
 9099 format (a1,1x,i3,19e10.5)
 9100 format (a1,1x,i3,5e10.5)

      end subroutine readcld

!***************************************************************************
      subroutine readaer(aerfile)
!***************************************************************************

      

! Purpose:  To read in IN_AER_RRTM, the file that contains input
!           aerosol properties.

! -------- Modules --------

      use rrsw_wvn, only : wavenum1, wavenum2

      implicit integer(i-n), real(a-h,o-z)
      
      character(len=*), intent(in) :: aerfile           ! input file name

! ------- Parameters -------
      parameter (mxlay = 203, jpbands = 29)
      parameter (ib1 = 16, ib2 = 29)
      parameter (mg = 16)
      parameter (mxstr = 16)
      parameter (mcmu = 32)

      real aerpar(3), ssa(jpbands), asym(jpbands), aod(mxlay),aod1(jpbands)
      real rlambda(jpbands), specfac(jpbands)
      real rnu0(16:29),rnu1(23:26)
      real f1(23:26),od0(23:26),od1(23:26)
      integer lay(mxlay),ivec(mxlay)

      common /control/ iaer, nstr, iout, istart, iend, icld, idelm, isccos
      common /profile/ nlayers,pavel(mxlay),tavel(mxlay), pz(0:mxlay),tz(0:mxlay),tbound
      common /swprop/  zenith, albedo, adjflux(jpbands)
      common /aerdat/  ssaaer(mxlay,jpbands), phase(mcmu,mxlay,jpbands), tauaer(mxlay,jpbands)

      character*1 ctest, cpercent

      data cpercent /'%'/

      data rnu0 /2903._8,3601._8,4310._8,4892._8,5623._8,6872._8, &
                 7872._8,10590._8,14420._8,18970._8,25015._8,30390._8, & 
                 43507._8,1412._8/
      data rnu1 /10530.7_8,14293.3_8,18678.0_8,24475.1_8/
 
      data f1  /0.9929_8,0.9883_8,0.978_8,0.9696_8/
      data od0 /0.1084_8,0.167_8,0.245_8,0.3611_8/
      data od1 /0.3144_8,0.4822_8,0.7013_8,1.0239_8/

      eps = 1.e-10_8
      irdaer = 12
      open(irdaer,file=aerfile,form='formatted')

      aod(:) = 0.0_8
      tauaer(:,ib1:ib2) = 0.0_8

! Read in number of different aerosol models option.
      read (irdaer, 9010) naer
!       if (naer .gt. 4) then
!          print *, 'NAER (= ', naer, ') IS GREATER THAN 4'
!          stop
!       endif
        
! For each aerosol read in optical properties and layer aerosol 
! optical depths.
      do ia = 1, naer
	 read (irdaer, 9011) nlay, iaod, issa, iasym, aerpar(1:3)

         if (iaod .eq. 0) then
! Set defaults to get standard Angstrom relation.
            if (aerpar(2) .lt. eps) aerpar(2) = 1._8

            do ib = ib1, ib2
  	       rlambda(ib) = 10000._8/rnu0(ib)
               specfac(ib) = (aerpar(2) + aerpar(3) * rlambda(ib)) / &
                            ((aerpar(2) + aerpar(3) - 1._8) + &
                              rlambda(ib)**aerpar(1))
            enddo
         endif

! For this aerosol, read in layers and optical depth information.
! Store a nonzero optical depth in aod to check for double specification.
         do il = 1, nlay
            read(irdaer, 9012) lay(il), aod1(ib1:ib2)
            if (aod(lay(il)) .lt. eps) then
               if (iaod .eq. 0) then
                  aod(lay(il)) = aod1(ib1)
                  do ib = ib1, ib2
                     tauaer(lay(il),ib) = aod(lay(il)) * specfac(ib)
                  enddo
               else
                  do ib = ib1, ib2
                     aod(lay(il)) = max(aod(lay(il)),aod1(ib))
                     tauaer(lay(il),ib) = aod1(ib)
                  enddo
               endif
            else
               print *,'LAYER ',lay(il),' HAS MORE THAN ONE AEROSOL TYPE'
               stop
            endif
         enddo

! Build vector of aerosol layer indices 

         do il=1,nlay
            ivec(il) = lay(il) 
         end do

! Correct bands 23 through 26 for sza effect (negligible for others)
         do ib=23,26
            if (iaod.eq.0) then
                od = sum(tauaer(ivec(1:nlay),ib))/zenith
                rnu = rnu0(ib) + &
                     (rnu1(ib)-rnu0(ib))*(od-od0(ib))/(od1(ib)-od0(ib))
               rlambda_new = 10000._8/rnu
               specfac_new = (aerpar(2)+aerpar(3)*rlambda_new) / &
                  ((aerpar(2)+aerpar(3)- 1.)+rlambda_new**aerpar(1))
               do il=1,nlay
                  tauaer(lay(il),ib) = tauaer(lay(il),ib) * &
                         specfac_new/specfac(ib)
               end do
            endif
         end do

! For this aerosol, read and store optical properties
         read (irdaer, 9013) ssa(ib1:ib2)

         do ib = ib1, ib2
            do il = 1, nlay
               if (issa .eq. 0) then 
                  ssaaer(lay(il),ib) = ssa(ib1)
               else
                  ssaaer(lay(il),ib) = ssa(ib)
               endif
            enddo
         enddo

         if (iasym .lt. 2) then
            read (irdaer, 9013) asym(ib1:ib2)

            do ib = ib1, ib2
               do il = 1, nlay
                  do istr = 1,  nstr
                     if (iasym .eq. 0) then 
                        phase(istr,lay(il),ib) = asym(ib1)**istr
                     elseif (iasym .eq. 1) then
                        phase(istr,lay(il),ib) = asym(ib)**istr
                     endif
                  enddo
               enddo
            enddo
         else
            do il = 1, nlay
               do istr = 1, nstr
                  read (irdaer, 9013) phase(istr,lay(il),ib1:ib2)
               enddo
            enddo
         endif

! End of naer loop
      enddo

 9000 continue
      close(irdaer)

 9010 format (3x, i2)
 9011 format (2x, i3, 4x, i1, 4x, i1, 4x, i1, 3f8.3)
 9012 format (2x, i3, 14f7.4)
 9013 format (14f5.3)

      end subroutine readaer

!**********************************************************************
      subroutine xsident(ird)
!**********************************************************************

! Purpose:  This subroutine identifies which cross-sections are to be used.

      implicit integer(i-n), real(a-h,o-z)

! ------- Parameters -------
      parameter (maxinpx=35)
      parameter (maxxsec=4)

      common /xsecctrl/ nxmol,ixindx(maxinpx)
                                                                         
!     nxmol     - number of cross-sections input by user
!     ixindx(i) - index of cross-section molecule corresponding to Ith
!                 cross-section specified by user
!                 = 0 -- not allowed in rrtm
!                 = 1 -- ccl4
!                 = 2 -- cfc11
!                 = 3 -- cfc12
!                 = 4 -- cfc22
!                                                                         
!     xsname=names, alias=aliases of the cross-section molecules          
!                                                                         
      character*10 xsname(maxinpx),alias(maxxsec,4),blank               

      data (alias(1,i),i=1,4)/ &
         'CCL4      ', 'CCL3F     ', 'CCL2F2    ', 'CHCLF2    '/ 
      data (alias(2,i),i=1,4)/ &
         ' ZZZZZZZZ ', 'CFCL3     ', 'CF2CL2    ', 'CHF2CL    '/         
      data (alias(3,i),i=1,4)/ &
         ' ZZZZZZZZ ', 'CFC11     ', 'CFC12     ', 'CFC22     '/         
      data (alias(4,i),i=1,4)/ &
         ' ZZZZZZZZ ', 'F11       ', 'F12       ', 'F22       '/        

      data blank / '          '/

      do i = 1, nxmol
         xsname(i) = blank
      enddo

! Read in the names of the molecules

      if (nxmol.gt.7) then                                               
         read (ird,'(7a10)') (xsname(i),i=1,7)                            
         read (ird,'(8a10)') (xsname(i),i=8,nxmol)                       
      else                                                                
         read (ird,'(7a10)') (xsname(i),i=1,nxmol)                       
      endif                                                               
                                                                         
!  Match the names read in against the names stored in alias           
!  and determine the index value.  
      ixmax = 4
      do i = 1, nxmol
!  Left-justify all inputed names.                                      
         call cljust (xsname(i),10)
         ixindx(i) = 0
         do j = 1, ixmax
            if ((xsname(i).eq.alias(1,j)) .or. &
                (xsname(i).eq.alias(2,j)) .or. &
                (xsname(i).eq.alias(3,j)) .or. &                           
                (xsname(i).eq.alias(4,j))) then                           
               ixindx(i) = j                                              
            endif                                                         
         enddo
      enddo   

      end subroutine xsident
      
end module read_input