program rrtmg_lw_test

      use read_input, only : readprof
      use rrtmg_lw_init, only: rrtmg_lw_ini
      use rrtmg_lw_rad, only: rrtmg_lw
      
      implicit none

! ------- Parameters -------      
      integer(kind=4), parameter :: mxlay = 203
      integer(kind=4), parameter :: mxmol = 38
      integer(kind=4), parameter :: nbndlw = 16
      integer(kind=4), parameter :: maxinpx = mxmol
      integer(kind=4), parameter :: maxxsec = 4
      integer(kind=4), parameter :: maxprod = mxlay*maxxsec
      integer(kind=8), parameter :: amd = 28.9660_8        ! Effective molecular weight of dry air (g/mol)
      integer(kind=8), parameter :: amw = 18.0160_8        ! Molecular weight of water vapor (g/mol)
      
! ------- Variables from Readprof -------      
      integer(kind=4) :: nlayers_out        ! total number of layers
      integer(kind=4) :: icld_out           ! clear/cloud flag
      integer(kind=4) :: imca               ! McICA on/off flag (1 = use McICA)
      integer(kind=4) :: iout_out           ! output option flag
      integer(kind=4) :: iaer_out           ! aerosol option flag
      integer(kind=4) :: idrv               ! Planck derivative option on/off flag 
                                                          ! (1 = provide upward flux adjustment
                                                          ! for change in surface temperature

      real(kind=8) :: pavel_out(mxlay)      ! layer pressures (mb) 
      real(kind=8) :: tavel_out(mxlay)      ! layer temperatures (K)
      real(kind=8) :: pz_out(0:mxlay)       ! level (interface) pressures (hPa, mb)
      real(kind=8) :: tz_out(0:mxlay)       ! level (interface) temperatures (K)
      real(kind=8) :: tbound_out            ! surface temperature (K)
      real(kind=8) :: dtbound_out           ! surface temperature change for idrv=1 (K)
      real(kind=8) :: coldry_out(mxlay)     ! dry air column density (mol/cm2)
      real(kind=8) :: wbrodl_out(mxlay)     ! broadening gas column density (mol/cm2)
      real(kind=8) :: wkl_out(mxmol,mxlay)  ! molecular amounts (mol/cm2)
      real(kind=8) :: wx_out(maxxsec,mxlay) ! cross-section amounts (mol/cm2)
      real(kind=8) :: pwvcm_out             ! precipitable water vapor (cm)
      real(kind=8) :: semiss_out(nbndlw)    ! lw surface emissivity

      integer(kind=4) :: inflag_out         ! cloud property option flag
      integer(kind=4) :: iceflag_out        ! ice cloud property flag
      integer(kind=4) :: liqflag_out        ! liquid cloud property flag

      real(kind=8) :: cldfrac_out(mxlay)    ! cloud fraction
      real(kind=8) :: tauc(nbndlw,mxlay)    ! in-cloud optical depth
!      real(kind=8) :: ssac(nbndlw,mxlay)   ! in-cloud single scattering albedo
                                                          !   for future expansion
!      real(kind=8) :: asmc(nbndlw,mxlay)   ! in-cloud asymmetry parameter
                                                          !   for future expansion
      real(kind=8) :: ciwp(mxlay)           ! in-cloud ice water path
      real(kind=8) :: clwp(mxlay)           ! in-cloud liquid water path
      real(kind=8) :: rel(mxlay)            ! cloud liquid particle effective radius (microns)
      real(kind=8) :: rei(mxlay)            ! cloud ice particle effective size (microns)
      real(kind=8) :: tauaer_out(mxlay,nbndlw)  ! aerosol optical depth
      
      character(len=64) :: filename, cldfile, aerfile

! ------- Variables for Computation -------  
      
      integer(kind=4) :: ncol, nlay, icld
      real(kind=8), allocatable :: play(:,:)
      real(kind=8), allocatable :: plev(:,:)
      real(kind=8), allocatable :: tlay(:,:)
      real(kind=8), allocatable :: tlev(:,:)
      real(kind=8), allocatable :: tsfc(:)
      real(kind=8), allocatable :: h2ovmr(:,:)
      real(kind=8), allocatable :: co2vmr(:,:)
      real(kind=8), allocatable :: o3vmr(:,:)
      real(kind=8), allocatable :: n2ovmr(:,:)
      real(kind=8), allocatable :: ch4vmr(:,:)
      real(kind=8), allocatable :: o2vmr(:,:)
      real(kind=8), allocatable :: cfc11vmr(:,:)
      real(kind=8), allocatable :: cfc12vmr(:,:)
      real(kind=8), allocatable :: cfc22vmr(:,:)
      real(kind=8), allocatable :: ccl4vmr(:,:)
      real(kind=8), allocatable :: emis(:,:)
      integer(kind=4) :: inflglw, iceflglw, liqflglw
      real(kind=8), allocatable :: cldfr(:,:)      
      real(kind=8), allocatable :: taucld(:,:,:)      
      real(kind=8), allocatable :: cicewp(:,:)      
      real(kind=8), allocatable :: cliqwp(:,:)      
      real(kind=8), allocatable :: reice(:,:)      
      real(kind=8), allocatable :: reliq(:,:)   
      real(kind=8), allocatable :: tauaer(:,:,:)     
      
      real(kind=8), allocatable :: uflx(:,:)
      real(kind=8), allocatable :: dflx(:,:)
      real(kind=8), allocatable :: hr(:,:)
      real(kind=8), allocatable :: uflxc(:,:)
      real(kind=8), allocatable :: dflxc(:,:)
      real(kind=8), allocatable :: hrc(:,:)
      real(kind=8), allocatable :: duflx_dt(:,:)
      real(kind=8), allocatable :: duflxc_dt(:,:)
      
      integer(kind=4) :: k
      
      call getarg(1, filename)
      call getarg(2, cldfile)
      call getarg(3, aerfile)
      

! ------- Main Program -------  
      call rrtmg_lw_ini(1.004e3)
      
      call readprof(nlayers_out, iout_out, imca, icld_out, &
          iaer_out, idrv, pavel_out, tavel_out, pz_out, tz_out, tbound_out, &
          semiss_out, dtbound_out, coldry_out, wkl_out, wbrodl_out, wx_out, &
          pwvcm_out, inflag_out, iceflag_out, liqflag_out, cldfrac_out, &
          tauc, ciwp, clwp, rei, rel, tauaer_out, filename, cldfile, aerfile)

      
      ! Allocate and prescribe values from readprof output
      ncol = 1
      nlay = nlayers_out
      icld = icld_out
      inflglw = inflag_out
      iceflglw = iceflag_out
      liqflglw = liqflag_out
      
      allocate(play    (ncol,nlay))
      allocate(plev    (ncol,nlay+1))
      allocate(tlay    (ncol,nlay))
      allocate(tlev    (ncol,nlay+1))
      allocate(tsfc    (ncol))
      allocate(h2ovmr  (ncol,nlay))
      allocate(o3vmr   (ncol,nlay))
      allocate(co2vmr  (ncol,nlay))
      allocate(ch4vmr  (ncol,nlay))
      allocate(n2ovmr  (ncol,nlay))
      allocate(o2vmr   (ncol,nlay))
      allocate(cfc11vmr(ncol,nlay))
      allocate(cfc12vmr(ncol,nlay))
      allocate(cfc22vmr(ncol,nlay))
      allocate(ccl4vmr (ncol,nlay))
      allocate(emis    (ncol,nbndlw))
      allocate(cldfr   (ncol,nlay))
      allocate(cliqwp  (ncol,nlay))
      allocate(cicewp  (ncol,nlay))
      allocate(reliq   (ncol,nlay))
      allocate(reice   (ncol,nlay))
      allocate(taucld  (nbndlw,ncol,nlay))
      allocate(tauaer  (ncol,nlay,nbndlw))
      
      allocate(uflx (ncol,nlay+1))
      allocate(dflx (ncol,nlay+1))
      allocate(hr   (ncol,nlay+1))
      allocate(uflxc(ncol,nlay+1))
      allocate(dflxc(ncol,nlay+1))
      allocate(hrc  (ncol,nlay+1))
      allocate(duflx_dt (ncol,nlay+1))
      allocate(duflxc_dt(ncol,nlay+1))
      
      ! p,t
      play(1,1:nlay) = pavel_out(1:nlay)
      plev(1,1:nlay+1) = pz_out(0:nlay)
      tlay(1,1:nlay) = tavel_out(1:nlay)
      tlev(1,1:nlay+1) = tz_out(0:nlay)
      tsfc(1) = tbound_out
      
      ! vmr
      h2ovmr(1,1:nlay) = wkl_out(1,1:nlay)/coldry_out(1:nlay)
      co2vmr(1,1:nlay) = wkl_out(2,1:nlay)/coldry_out(1:nlay) 
      o3vmr (1,1:nlay) = wkl_out(3,1:nlay)/coldry_out(1:nlay) 
      n2ovmr(1,1:nlay) = wkl_out(4,1:nlay)/coldry_out(1:nlay) 
      ch4vmr(1,1:nlay) = wkl_out(6,1:nlay)/coldry_out(1:nlay)
      o2vmr (1,1:nlay) = wkl_out(7,1:nlay)/coldry_out(1:nlay)
      
      ccl4vmr (1,1:nlay) = wx_out(1,1:nlay)/coldry_out(1:nlay)*1.e+20
      cfc11vmr(1,1:nlay) = wx_out(2,1:nlay)/coldry_out(1:nlay)*1.e+20
      cfc12vmr(1,1:nlay) = wx_out(3,1:nlay)/coldry_out(1:nlay)*1.e+20
      cfc22vmr(1,1:nlay) = wx_out(4,1:nlay)/coldry_out(1:nlay)*1.e+20
      emis  (1,1:nbndlw) = semiss_out(1:nbndlw)
      
      !cloud, aerosol ...
      cldfr  (1,1:nlay)  = cldfrac_out(1:nlay)
      cicewp (1,1:nlay)  = ciwp (1:nlay)
      cliqwp (1,1:nlay)  = clwp (1:nlay)
      reice  (1,1:nlay)  = rei  (1:nlay)
      reliq  (1,1:nlay)  = rel  (1:nlay)
      do k = 1, nbndlw
          taucld(k,1,1:nlay) = tauc(k, 1:nlay)
          tauaer(1,1:nlay,k) = tauaer_out(1:nlay,k)
      end do
      
      ! This is for No_MCICA only...
      call rrtmg_lw(ncol, nlay, icld, idrv, play,plev,tlay,tlev,tsfc,&
                    h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr, &
                    cfc11vmr, cfc12vmr, cfc22vmr, ccl4vmr, emis, &
                    inflglw, iceflglw, liqflglw, cldfr, &
                    taucld, cicewp, cliqwp, reice, reliq, tauaer, &
                    uflx, dflx, hr(:,1:nlay), uflxc, dflxc, hrc(:,1:nlay), duflx_dt,duflxc_dt )
      hr (1,nlay+1) = 0.0_8
      hrc(1,nlay+1) = 0.0_8
                    
      ! Output
      
      write(*,9900)
      write(*,9901)

      do k = nlay+1, 1, -1
            if (plev(1,k) .lt. 1.e-2_8) then
                write(*,9952) k-1,plev(1,k),uflx(1,k),dflx(1,k),uflx(1,k)-dflx(1,k),hr(1,k)
            elseif (plev(1,k) .lt. 1.e-1) then
                write(*,9953) k-1,plev(1,k),uflx(1,k),dflx(1,k),uflx(1,k)-dflx(1,k),hr(1,k)
            elseif (plev(1,k) .lt. 1.) then
                write(*,9954) k-1,plev(1,k),uflx(1,k),dflx(1,k),uflx(1,k)-dflx(1,k),hr(1,k)
            elseif (plev(1,k) .lt. 10.) then
                write(*,9955) k-1,plev(1,k),uflx(1,k),dflx(1,k),uflx(1,k)-dflx(1,k),hr(1,k)
            elseif (plev(1,k) .lt. 100.) then
                write(*,9956) k-1,plev(1,k),uflx(1,k),dflx(1,k),uflx(1,k)-dflx(1,k),hr(1,k)
            elseif (plev(1,k) .lt. 1000.) then
                write(*,9957) k-1,plev(1,k),uflx(1,k),dflx(1,k),uflx(1,k)-dflx(1,k),hr(1,k)
            else
                write(*,9958) k-1,plev(1,k),uflx(1,k),dflx(1,k),uflx(1,k)-dflx(1,k),hr(1,k)
            endif
      enddo
      
 9952 format(1x,i3,9x,f7.6,3x,f8.4,6x,f8.4,6x,f12.7,10x,f9.5)
 9953 format(1x,i3,9x,f6.5,4x,f8.4,6x,f8.4,6x,f12.7,10x,f9.5)
 9954 format(1x,i3,8x,f6.4,5x,f8.4,6x,f8.4,6x,f12.7,10x,f9.5)
 9955 format(1x,i3,7x,f6.3,6x,f8.4,6x,f8.4,6x,f12.7,10x,f9.5)
 9956 format(1x,i3,6x,f6.2,7x,f8.4,6x,f8.4,6x,f12.7,10x,f9.5)
 9957 format(1x,i3,5x,f6.1,8x,f8.4,6x,f8.4,6x,f12.7,10x,f9.5)
 9958 format(1x,i3,5x,f6.1,8x,f8.4,6x,f8.4,6x,f12.7,10x,f9.5)
 9900 format(1x,'LEVEL    PRESSURE   UPWARD FLUX   DOWNWARD FLUX    NET FLUX       HEATING RATE')
 9901 format(1x,'            mb          W/m2          W/m2           W/m2          degree/day')





end program rrtmg_lw_test