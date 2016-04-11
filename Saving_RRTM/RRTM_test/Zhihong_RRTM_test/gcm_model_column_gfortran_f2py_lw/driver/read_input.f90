module read_input

! -------- Modules --------

      implicit none
      public :: readprof, readcld, readaer

      contains

!************************************************************************
      subroutine readprof(nlayers_out, iout_out, imca, icld_out, &
          iaer_out, idrv, pavel_out, tavel_out, pz_out, tz_out, tbound_out, &
          semiss_out, dtbound_out, coldry_out, wkl_out, wbrodl_out, wx_out, &
          pwvcm_out, inflag_out, iceflag_out, liqflag_out, cldfrac_out, &
          tauc, ciwp, clwp, rei, rel, tauaer_out, filename, cldfile, aerfile)
!************************************************************************

! --------- Modules ----------

      use rrlw_con, only: pi, grav, planck, boltz, clight, avogad, alosmt, &
                          gascon, radcn1, radcn2 

! Note: COMMON blocks are left in this routine for array passing with 
!       rrtatm.f for reading input profiles in single column mode.
!       Scalars and arrays in subroutine call are renamed to avoid conflict
!       with COMMON blocks. Variable mxlay is the maximum possible number of
!       layers, which is used to dimension the input arrays, while nlayers
!       is the actual number of input layers. 
!
! Purpose: Read in atmospheric profile.

      implicit integer(i-n), real (a-h,o-z)

! ------- Parameters -------
      parameter (mxlay = 203)
      parameter (mxmol = 38)
      parameter (nbndlw = 16)
      parameter (maxinpx = mxmol)
      parameter (maxxsec = 4)
      parameter (maxprod = mxlay*maxxsec)
      parameter (amd = 28.9660_8)         ! Effective molecular weight of dry air (g/mol)
      parameter (amw = 18.0160_8)         ! Molecular weight of water vapor (g/mol)

      dimension altz(0:mxlay),ixtrans(14),semis(nbndlw)

      common /consts/   pic,planckc,boltzc,clightc,avogadc,alosmtc,gasconc, &
                        radcn1c,radcn2c
      common /control/  numangs, iout, istart, iend, icld, iaer
      common /profile/  nlayers,pavel(mxlay),tavel(mxlay),pz(0:mxlay),tz(0:mxlay)
      common /surface/  tbound,ireflect,semiss(nbndlw)
      common /species/  coldry(mxlay),wkl(mxmol,mxlay),wbrodl(mxlay),colmol(mxlay),nmol
      common /ifil/     ird,ipr,ipu,idum(15)
      common /xsecctrl/ nxmol,ixindx(maxinpx)
      common /xsec/     wx(maxxsec,mxlay)
      common /pathx/    ixmax,nxmol0,ixindx0(maxinpx),wx0(maxinpx,mxlay)    
      common /xrrtatm/  ixsect

      common /cloudin/   inflag,clddat1(mxlay),clddat2(mxlay), &
                         iceflag,liqflag,clddat3(mxlay),clddat4(mxlay)
      common /clouddat/  ncbands,cldfrac(mxlay),taucloud(mxlay,nbndlw)
      common /aerdat/    tauaer(mxlay,nbndlw)

      character*80 form1(0:1),form2(0:1),form3(0:1)
      character*1 ctest, cdollar, cprcnt,cdum

! Dimensions for transfer to rrtmg
      ! integer(kind=4), intent(in) :: ird_in              ! input file unit
      character(len=*), intent(in) :: filename           ! input file name
      character(len=*), intent(in), optional :: cldfile  ! input file name
      character(len=*), intent(in), optional :: aerfile  ! input file name
      integer(kind=4), intent(out) :: nlayers_out        ! total number of layers
      integer(kind=4), intent(out) :: icld_out           ! clear/cloud flag
      integer(kind=4), intent(out) :: imca               ! McICA on/off flag (1 = use McICA)
      integer(kind=4), intent(out) :: iout_out           ! output option flag
      integer(kind=4), intent(out) :: iaer_out           ! aerosol option flag
      integer(kind=4), intent(out) :: idrv               ! Planck derivative option on/off flag 
                                                          ! (1 = provide upward flux adjustment
                                                          ! for change in surface temperature

      real(kind=8), intent(out) :: pavel_out(mxlay)      ! layer pressures (mb) 
      real(kind=8), intent(out) :: tavel_out(mxlay)      ! layer temperatures (K)
      real(kind=8), intent(out) :: pz_out(0:mxlay)       ! level (interface) pressures (hPa, mb)
      real(kind=8), intent(out) :: tz_out(0:mxlay)       ! level (interface) temperatures (K)
      real(kind=8), intent(out) :: tbound_out            ! surface temperature (K)
      real(kind=8), intent(out) :: dtbound_out           ! surface temperature change for idrv=1 (K)
      real(kind=8), intent(out) :: coldry_out(mxlay)     ! dry air column density (mol/cm2)
      real(kind=8), intent(out) :: wbrodl_out(mxlay)     ! broadening gas column density (mol/cm2)
      real(kind=8), intent(out) :: wkl_out(mxmol,mxlay)  ! molecular amounts (mol/cm2)
      real(kind=8), intent(out) :: wx_out(maxxsec,mxlay) ! cross-section amounts (mol/cm2)
      real(kind=8), intent(out) :: pwvcm_out             ! precipitable water vapor (cm)
      real(kind=8), intent(out) :: semiss_out(nbndlw)    ! lw surface emissivity

      integer(kind=4), intent(out) :: inflag_out         ! cloud property option flag
      integer(kind=4), intent(out) :: iceflag_out        ! ice cloud property flag
      integer(kind=4), intent(out) :: liqflag_out        ! liquid cloud property flag

      real(kind=8), intent(out) :: cldfrac_out(mxlay)    ! cloud fraction
      real(kind=8), intent(out) :: tauc(nbndlw,mxlay)    ! in-cloud optical depth
!      real(kind=8), intent(out) :: ssac(nbndlw,mxlay)   ! in-cloud single scattering albedo
                                                          !   for future expansion
!      real(kind=8), intent(out) :: asmc(nbndlw,mxlay)   ! in-cloud asymmetry parameter
                                                          !   for future expansion
      real(kind=8), intent(out) :: ciwp(mxlay)           ! in-cloud ice water path
      real(kind=8), intent(out) :: clwp(mxlay)           ! in-cloud liquid water path
      real(kind=8), intent(out) :: rel(mxlay)            ! cloud liquid particle effective radius (microns)
      real(kind=8), intent(out) :: rei(mxlay)            ! cloud ice particle effective size (microns)
      real(kind=8), intent(out) :: tauaer_out(mxlay,nbndlw)  ! aerosol optical depth
!      real(kind=8), intent(out) :: ssaaer_out(mxlay,nbndlw)  ! aerosol single scattering albedo
                                                               !   for future expansion
!      real(kind=8), intent(out) :: asmaer_out(mxlay,nbndlw)  ! aerosol asymmetry parameter
                                                               !   for future expansion

! Local
      real(kind=8) :: fice(mxlay)                        ! cloud ice fraction

      real(kind=8) :: amttl                              ! moist air vertical sum (molecular amount)
      real(kind=8) :: wvttl                              ! water vapor vertical sum (molecular amount)
      real(kind=8) :: summol                             ! sum over non-water molecules
      real(kind=8) :: wvsh                               ! water vapor vertical total specific humitidy
      real(kind=8) :: pwvcm                              ! precipitable water vapor (cm)
      real(kind=8) :: dtbound                            ! change in surface temperature for idrv=1 (K)


!

! Initializations

      data cdollar /'$'/
      data cprcnt /'%'/
      data ixtrans /0,0,0,1,2,3,0,0,0,0,0,4,0,0/

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

      do l = 1, mxlay
         do ix = 1, maxxsec
           wx(ix,l) = 0.0_8
         enddo
      enddo

! Top of read input loop
 1000 continue
      read (ird,9010,end=8800) ctest
      if (ctest .eq. cprcnt) goto 8900 
      if (ctest .ne. cdollar) goto 1000

      read (ird,9011) iaer, iatm, ixsect, numangs, iout, idrv, imca, icld

!  If numangs set to -1, reset to default rt code for
!  backwards compatibility with original rrtm
      if (numangs .eq. -1) numangs = 0

!  If clouds are present, read in appropriate input file, IN_CLD_RRTM.
      if (icld .ge. 1) call readcld(cldfile)

!  If aerosols are present, read in appropriate input file, IN_AER_RRTM.
      if (iaer .eq. 10) call readaer(aerfile)

!  Read in surface information.
      read (ird,9012) tbound,iemiss,ireflect,semis(1:nbndlw)

!  Read in change in surface temperature for upward flux derivative adjustment 
      dtbound = 0.0_8
      if (idrv .eq. 1) then
         read (ird,9012) dtbound
      endif

      do iband = 1, nbndlw
         semiss(iband) = 1.0_8
         if (iemiss .eq. 1 .and. semis(1) .ne. 0._8) then
            semiss(iband) = semis(1)
         elseif (iemiss .eq. 2) then
            if (semis(iband) .ne. 0._8) then
               semiss(iband) = semis(iband)
            endif
         endif
      enddo

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
               if (nxmol0 .gt. 7) read (ird,form3(iformx)) (wx0(m,l),m=8,nxmol0)
            enddo
         endif
      else
         ipu = 7
         ipr = 66
         open(unit=ipr,file='TAPE6',status='unknown')
         call rrtatm
         if (ixsect .eq. 1) then
            do mx = 1, nxmol0
               ixindx(mx) = ixtrans(ixindx0(mx))
            enddo
         endif
      endif
      if (tbound .lt. 0) tbound = tz(0)

!  Test for mixing ratio input.
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
      amttl = 0.0_8
      wvttl = 0.0_8
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
         amttl = amttl + coldry(l)+wkl(1,l)
         wvttl = wvttl + wkl(1,l)
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

!  Calculate total precipitable water 
      wvsh = (amw * wvttl) / (amd * amttl)
      pwvcm = wvsh * (1.e3_8 * pz(0)) / (1.e2_8 * grav)

! Pass output arrays to new variables for transfer to rrtmg through subroutine call.
      nlayers_out = nlayers
      iout_out = iout
      icld_out = icld
      iaer_out = iaer
      tbound_out = tbound
      dtbound_out = dtbound
      pwvcm_out = pwvcm
      inflag_out = inflag
      iceflag_out = iceflag
      liqflag_out = liqflag

      pz_out(0) = pz(0)
      tz_out(0) = tz(0)
      do l = 1, nlayers
         pavel_out(l) = pavel(l)
         tavel_out(l) = tavel(l)
         pz_out(l) = pz(l)
         tz_out(l) = tz(l)
         coldry_out(l) = coldry(l)
         wbrodl_out(l) = wbrodl(l)
         cldfrac_out(l) = cldfrac(l)
         do m = 1, mxmol
            wkl_out(m,l) = wkl(m,l)
         enddo  
         do ix = 1, maxxsec
            wx_out(ix,l) = wx(ix,l)
         enddo  
      enddo
      do n = 1, nbndlw
         semiss_out(n) = semiss(n)
      enddo

      do l = 1, nlayers
         if (inflag.eq.0) then
            do n = 1, nbndlw
               tauc(n,l) = clddat1(l)
!               ssac(n,l) = 1._8
!               asmc(n,l) = 1._8
            enddo
            ciwp(l) = 0._8
            clwp(l) = 0._8
            fice(l) = 0._8
            rei(l) = 0._8
            rel(l) = 0._8
         else
            do n = 1, nbndlw
               tauc(n,l) = 0._8
!               ssac(n,l) = 1._8
!               asmc(n,l) = 1._8
            enddo
            cwp = clddat1(l)
            fice(l) = clddat2(l)
            ciwp(l) = cwp * fice(l)
            clwp(l) = cwp * (1._8 - fice(l))
            rei(l) = clddat3(l)
            rel(l) = clddat4(l)
         endif 
      enddo

      do l = 1, nlayers
         do n = 1, nbndlw
            tauaer_out(l,n) = tauaer(l,n)
!            ssaaer_out(l,n) = 1._8
!            asmaer_out(l,n) = 0._8
         enddo
      enddo

      goto 9000

 8800 continue
 8900 if (ctest.eq.'%') stop 'END OF INPUT FILE'
 9000 continue

 9010 format (a1)
 9011 format (18x,i2,29x,i1,19x,i1,13x,i2,2x,i3,1x,i1,1x,i1,i1)
 9012 format (e10.3,1x,i1,2x,i1,16e5.3)
 9013 format (1x,i1,i3,i5)                                     
 9300 format (i5)
 9301 format (1x,i1)

      end subroutine readprof

      
!*************************************************************************
      subroutine readcld(cldfile)
!*************************************************************************

! --------- Modules ----------

      
! Purpose:  To read in IN_CLD_RRTM, the file that contains input 
!           cloud properties.
      implicit integer(i-n), real (a-h,o-z)

      character(len=*), intent(in) :: cldfile           ! input file name
      
! ------- Parameters -------
      parameter (mxlay=203)
      parameter (nbndlw = 16)

      common /profile/   nlayers,pavel(mxlay),tavel(mxlay),pz(0:mxlay),tz(0:mxlay)
      common /cloudin/   inflag,clddat1(mxlay),clddat2(mxlay), &
                         iceflag,liqflag,clddat3(mxlay),clddat4(mxlay)
      common /clouddat/  ncbands,cldfrac(mxlay),taucloud(mxlay,nbndlw)

      character*1 ctest, cpercent

      data cpercent /'%'/
      irdcld = 11

      open(irdcld,file=cldfile,form='formatted')

! Read in cloud input option.  
      read(irdcld,9050) inflag, iceflag, liqflag

      do lay = 1, nlayers
         cldfrac(lay) = 0._8
      enddo

! Top of read input loop
 1000 continue

!  For INFLAG = 0 or 1, for each cloudy layer only LAY, FRAC, and
!  DAT1 are pertinent.  If CTEST = '%', then there are no more 
!  cloudy layers to process.
      read (irdcld,9100,end=9000) ctest,lay,frac,dat1,dat2,dat3,dat4
      if (ctest .eq. cpercent) goto 9000
      cldfrac(lay) = frac
      clddat1(lay) = dat1
      clddat2(lay) = dat2
      clddat3(lay) = dat3
      clddat4(lay) = dat4
      goto 1000

 9000 continue
      close(irdcld)

 9050 format (3x,i2,4x,i1,4x,i1)
 9100 format (a1,1x,i3,5e10.5)

      end subroutine readcld

!***************************************************************************
      subroutine readaer(aerfile)
!***************************************************************************


! Purpose:  To read in IN_AER_RRTM, the file that contains input
!           aerosol properties.

! -------- Modules --------

      use rrlw_wvn, only : wavenum1, wavenum2
      
      implicit integer(i-n), real(a-h,o-z)
      
      character(len=*), intent(in) :: aerfile           ! input file name


! ------- Parameters -------
      parameter (mxlay = 203)
      parameter (nbndlw  = 16)
!      parameter (mg = 16)
!      parameter (mxstr = 16)
!      parameter (mcmu = 32)

      real aod(mxlay),aod1(nbndlw)
      integer lay(mxlay),ivec(mxlay)

      common /control/  numangs, iout, istart, iend, icld, iaer
      common /profile/  nlayers,pavel(mxlay),tavel(mxlay),pz(0:mxlay),tz(0:mxlay)

      common /aerdat/  tauaer(mxlay,nbndlw)

      character*1 ctest, cpercent

      data cpercent /'%'/

      eps = 1.e-10_8
      irdaer = 12
      open(irdaer,file=aerfile,form='formatted')

      aod(:) = 0.0_8
      tauaer(:,:) = 0.0_8

! Read in number of different aerosol models option.
      read (irdaer, 9010) naer
!       if (naer .gt. 4) then
!          print *, 'NAER (= ', naer, ') IS GREATER THAN 4'
!          stop
!       endif
        
! For each aerosol read in optical properties and layer aerosol 
! optical depths.
      do ia = 1, naer
	 read (irdaer, 9011) nlay, iaod

! Input restricted to direct input of aerosol optical depths
         iaod = 1

! For this aerosol, read in layers and optical depth information.
! Store a nonzero optical depth in aod to check for double specification.
         do il = 1, nlay
            read(irdaer, 9012) lay(il), (aod1(ib), ib = 1,nbndlw)
            if (aod(lay(il)) .lt. eps) then
               if (iaod .eq. 1) then
                  do ib = 1, nbndlw
                     aod(lay(il)) = max(aod(lay(il)),aod1(ib))
                     tauaer(lay(il),ib) = aod1(ib)
                  enddo
               endif
            else
               print *,'LAYER ',lay(il),' HAS MORE THAN ONE AEROSOL TYPE'
               stop
            endif
         enddo

! End of naer loop
      enddo

 9000 continue
      close(irdaer)

 9010 format (3x, i2)
 9011 format (2x, i3, 4x, i1)
 9012 format (2x, i3, 16f7.4)

      end subroutine readaer
      
      
!**********************************************************************
      subroutine xsident(ird)
!**********************************************************************

! Purpose:  This subroutine identifies which cross-sections are to be used.

      implicit integer(i-n), real (a-h,o-z)

! ------- Parameters -------
      parameter (maxinpx=38)
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