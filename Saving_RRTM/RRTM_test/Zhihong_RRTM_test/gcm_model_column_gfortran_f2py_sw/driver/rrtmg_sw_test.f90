program rrtmg_sw_test

       
      use read_input, only : readprof
      use rrtmg_sw_init, only: rrtmg_sw_ini
      use rrtmg_sw_rad, only: rrtmg_sw
      
      implicit none

! ------- Parameters from parrrsw.f90 -------     
      integer(kind=4), parameter :: mxlay  = 203    !jplay, klev
      integer(kind=4), parameter :: nbndsw = 14     !jpsw, ksw
      integer(kind=4), parameter :: naerec  = 6     !jpaer
      integer(kind=4), parameter :: mxmol  = 38
      integer(kind=4), parameter :: nstr   = 2
      integer(kind=4), parameter :: nmol   = 7
! Use for 112 g-point model   
      integer(kind=4), parameter :: ngptsw = 112    !jpgpt
! Use for 224 g-point model   
!      integer(kind=4), parameter :: ngptsw = 224   !jpgpt

! may need to rename these - from v2.6
      integer(kind=4), parameter :: jpbands   = 29
      integer(kind=4), parameter :: jpb1     = 16   !istart
      integer(kind=4), parameter :: jpb2     = 29   !iend

! Source function solar constant
      real(kind=8), parameter :: rrsw_scon = 1.36822e+03     ! W/m2 
      
! ------- Variables from Readprof -------      
      integer(kind=4) :: nlayers_out           ! total number of layers
      integer(kind=4) :: imca                  ! McICA on/off flag (1 = use McICA)
      integer(kind=4) :: icld_out              ! clear/cloud/overlap flag
      integer(kind=4) :: iout_out              ! output option flag
      integer(kind=4) :: iaer_out              ! aerosol flag
      integer(kind=4) :: isccos_out            ! aerosol flag
      integer(kind=4) :: idelm_out             ! aerosol flag

      real(kind=8) :: pavel_out(mxlay)         ! layer pressures (mb) 
      real(kind=8) :: tavel_out(mxlay)         ! layer temperatures (K)
      real(kind=8) :: pz_out(0:mxlay)          ! level (interface) pressures (hPa, mb)
      real(kind=8) :: tz_out(0:mxlay)          ! level (interface) temperatures (K)
      real(kind=8) :: tbound_out               ! surface temperature (K)
      real(kind=8) :: pdp(mxlay)               ! layer pressure thickness (hPa, mb)
      real(kind=8) :: coldry_out(mxlay)        ! dry air molecular amount
      real(kind=8) :: wkl_out(mxmol,mxlay)     ! molecular amounts (mol/cm-2)
      real(kind=8) :: semiss_out(jpbands)       ! surface emissivity
      real(kind=8) :: zenith_out               ! cos solar zenith angle
      real(kind=8) :: adjflux_out(jpbands)      ! adjustment for current Earth/Sun distance

      integer(kind=4) :: inflag_out              ! cloud property option flag
      integer(kind=4) :: iceflag_out             ! ice cloud property flag
      integer(kind=4) :: liqflag_out             ! liquid cloud property flag

      real(kind=8) :: cldfrac_out(mxlay)         ! cloud fraction
      real(kind=8) :: tauc(nbndsw,mxlay)         ! in-cloud optical depth (non-delta scaled)
      real(kind=8) :: ssac(nbndsw,mxlay)         ! in-cloud single scattering albedo (non-delta scaled)
      real(kind=8) :: asmc(nbndsw,mxlay)         ! in-cloud asymmetry parameter (non-delta scaled)
      real(kind=8) :: fsfc(nbndsw,mxlay)         ! in-cloud forward scattering fraction (non-delta scaled)
      real(kind=8) :: ciwp(mxlay)                ! in-cloud ice water path
      real(kind=8) :: clwp(mxlay)                ! in-cloud liquid water path
      real(kind=8) :: rei(mxlay)                 ! cloud ice particle size
      real(kind=8) :: rel(mxlay)                 ! cloud liquid particle size
      real(kind=8) :: tauaer_out(mxlay,jpbands)   ! aerosol optical depth
      real(kind=8) :: ssaaer_out(mxlay,jpbands)   ! aerosol single scattering albedo
      real(kind=8) :: asmaer_out(mxlay,jpbands)   ! aerosol asymmetry parameter
                                                                 !   first momemnt of input phase function
      ! Added by ZTAN
      integer(kind=4) :: dyofyr_out
      real(kind=8) :: adjes_out
      character(len=64) :: filename, cldfile, aerfile
      
      
! ------- Variables for Computation -------  
      
      integer(kind=4) :: ncol, nlay, icld, iaer
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

      real(kind=8), allocatable :: asdir(:) 
      real(kind=8), allocatable :: aldir(:)
      real(kind=8), allocatable :: asdif(:)
      real(kind=8), allocatable :: aldif(:)
      
      integer(kind=4) :: dyofyr
      real(kind=8) :: adjes, scon
      real(kind=8), allocatable :: coszen(:)
      
      integer(kind=4) :: inflgsw, iceflgsw, liqflgsw
      real(kind=8), allocatable :: cldfr(:,:)      
      real(kind=8), allocatable :: taucld(:,:,:)    
      real(kind=8), allocatable :: ssacld(:,:,:)   
      real(kind=8), allocatable :: asmcld(:,:,:)   
      real(kind=8), allocatable :: fsfcld(:,:,:)
      real(kind=8), allocatable :: cicewp(:,:)      
      real(kind=8), allocatable :: cliqwp(:,:)      
      real(kind=8), allocatable :: reice(:,:)      
      real(kind=8), allocatable :: reliq(:,:)   
      real(kind=8), allocatable :: tauaer(:,:,:)  
      real(kind=8), allocatable :: ssaaer(:,:,:) 
      real(kind=8), allocatable :: asmaer(:,:,:) 
      real(kind=8), allocatable :: ecaer(:,:,:)    
      
      real(kind=8), allocatable :: uflx(:,:)
      real(kind=8), allocatable :: dflx(:,:)
      real(kind=8), allocatable :: hr(:,:)
      real(kind=8), allocatable :: uflxc(:,:)
      real(kind=8), allocatable :: dflxc(:,:)
      real(kind=8), allocatable :: hrc(:,:)
      
      real(kind=8), allocatable :: dirdflux(:,:)
      real(kind=8), allocatable :: difdflux(:,:)
      
      integer(kind=4) :: k
      integer(kind=4) :: indform, idelm
      character page
      character*50 outform(7)
      
      data outform &
         /'(1x,i3,3x,f7.6,4x,4(f10.4,4x),f11.6,4x,f10.5)',&
          '(1x,i3,4x,f6.5,4x,4(f10.4,4x),f11.6,4x,f10.5)',&
          '(1x,i3,4x,f6.4,4x,4(f10.4,4x),f11.6,4x,f10.5)',&
          '(1x,i3,4x,f6.3,4x,4(f10.4,4x),f11.6,4x,f10.5)',&
          '(1x,i3,4x,f6.2,4x,4(f10.4,4x),f11.6,4x,f10.5)',&
          '(1x,i3,4x,f6.1,4x,4(f10.4,4x),f11.6,4x,f10.5)',&
          '(1x,i3,4x,f6.1,4x,4(f10.4,4x),f11.6,4x,f10.5)'/ 
          
      page = char(12)
      
      call getarg(1, filename)
      call getarg(2, cldfile)
      call getarg(3, aerfile)
      

! ------- Main Program -------  

      call rrtmg_sw_ini(1.004e3)
      
      call readprof(nlayers_out, iout_out, imca, icld_out, &
           iaer_out, isccos_out, idelm_out, pdp, &
           pavel_out, tavel_out, pz_out, tz_out, tbound_out, semiss_out, &
           zenith_out, adjflux_out, dyofyr_out, adjes_out, &
           coldry_out, wkl_out, inflag_out, iceflag_out,liqflag_out, &
           cldfrac_out, tauc, ssac, asmc, fsfc, ciwp, clwp, rei, rel, &
           tauaer_out, ssaaer_out, asmaer_out, filename, cldfile, aerfile)
      ! write(*,*), 'idelm from input: ', idelm_out
      
      
      ! Allocate and prescribe values from readprof output
      ncol = 1
      nlay = nlayers_out
      icld = icld_out
      iaer = iaer_out
      inflgsw = inflag_out
      iceflgsw = iceflag_out
      liqflgsw = liqflag_out
      dyofyr   = dyofyr_out
      adjes    = adjes_out
      scon     = rrsw_scon 
      idelm    = idelm_out
      
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
      
      allocate(asdir   (ncol)) 
      allocate(aldir   (ncol))
      allocate(asdif   (ncol))
      allocate(aldif   (ncol))
      allocate(coszen  (ncol))
      
      allocate(cldfr   (ncol,nlay))
      allocate(taucld  (nbndsw,ncol,nlay))
      allocate(ssacld  (nbndsw,ncol,nlay))
      allocate(asmcld  (nbndsw,ncol,nlay))
      allocate(fsfcld  (nbndsw,ncol,nlay))
      
      allocate(cicewp  (ncol,nlay))
      allocate(cliqwp  (ncol,nlay))
      allocate(reliq   (ncol,nlay))
      allocate(reice   (ncol,nlay))
      
      allocate(tauaer  (ncol,nlay,nbndsw))
      allocate(ssaaer  (ncol,nlay,nbndsw))
      allocate(asmaer  (ncol,nlay,nbndsw))
      allocate(ecaer   (ncol,nlay,naerec))
      
      allocate(uflx (ncol,nlay+1))
      allocate(dflx (ncol,nlay+1))
      allocate(hr   (ncol,nlay+1))
      allocate(uflxc(ncol,nlay+1))
      allocate(dflxc(ncol,nlay+1))
      allocate(hrc  (ncol,nlay+1))
      
      allocate(dirdflux (ncol,nlay+1))
      allocate(difdflux (ncol,nlay+1))
      
      
      ! ___ COMPLETED ABOVE ___ !    
      
      
      
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
      
      ! asdir, aldir, asdif, aldif <- semiss
      asdir(1) = 1._8 - semiss_out(jpb1)
      aldir(1) = 1._8 - semiss_out(jpb1)
      asdif(1) = 1._8 - semiss_out(jpb1)
      aldif(1) = 1._8 - semiss_out(jpb1)
      ! Note: albedo by band is not implemented !!
      
      ! coszen
      coszen(1) = zenith_out
      ! write(*,*)  zenith_out
      
      !cloud, aerosol ...
      cldfr  (1,1:nlay)  = cldfrac_out(1:nlay)
      cicewp (1,1:nlay)  = ciwp (1:nlay)
      cliqwp (1,1:nlay)  = clwp (1:nlay)
      reice  (1,1:nlay)  = rei  (1:nlay)
      reliq  (1,1:nlay)  = rel  (1:nlay)
      
      ! cloud and aerosol 
      do k = jpb1,jpb2
          ! write(*,*),k, k-jpb1+1
          taucld(k-jpb1+1,1,1:nlay) = tauc(k, 1:nlay)
          ssacld(k-jpb1+1,1,1:nlay) = ssac(k, 1:nlay)
          asmcld(k-jpb1+1,1,1:nlay) = asmc(k, 1:nlay)
          fsfcld(k-jpb1+1,1,1:nlay) = fsfc(k, 1:nlay)
          
          tauaer(1,1:nlay,k-jpb1+1) = tauaer_out(1:nlay,k)
          ssaaer(1,1:nlay,k-jpb1+1) = ssaaer_out(1:nlay,k)
          asmaer(1,1:nlay,k-jpb1+1) = asmaer_out(1:nlay,k)
      end do
      
      !write(*,*), iaer
      !write(*,*), tauaer_out(1:12,jpb1:jpb2)
      !write(*,*), ssaaer_out(1:12,jpb1:jpb2)
      !write(*,*), asmaer_out(1:12,jpb1:jpb2)
      
      ! ecaer
      ecaer(1,1:nlay,1:naerec) = 1.0e-15_8
      
      ! This is for No_MCICA only...
      call rrtmg_sw(ncol, nlay, icld, iaer, idelm, play,plev,tlay,tlev,tsfc,&
                    h2ovmr, o3vmr, co2vmr, ch4vmr, n2ovmr, o2vmr, &
                    asdir   ,asdif   ,aldir   ,aldif   , &  ! Added for SW
                    coszen  ,adjes   ,dyofyr  ,scon    , &  ! Added for SW
                    inflgsw, iceflgsw, liqflgsw, cldfr, &
                    taucld, ssacld, asmcld, fsfcld, & ! ssacld, asmcld, fsfcld: Added for SW
                    cicewp, cliqwp, reice, reliq, & 
                    tauaer, ssaaer, asmaer, ecaer, &  ! ssaaer, asmaer, ecaer: Added for SW
                    uflx, dflx, hr(:,1:nlay), uflxc, dflxc, hrc(:,1:nlay), &
                    dirdflux, difdflux)               ! dirdflux, difdflux: Added for SW
                    
      ! deleted: cfc11vmr, cfc12vmr, cfc22vmr, ccl4vmr, emis, duflx_dt,duflxc_dt 
      
      hr (1,nlay+1) = 0.0_8
      hrc(1,nlay+1) = 0.0_8
                  
      ! Output
      
      
      ! ___ COMPLETED BELOW ___ !
      ! ___ OUTPUT of isccos, idelm are not included yet ___ !
      write(*,9900)
      write(*,9901)

      do k = nlay+1, 1, -1
         if (plev(1,k) .lt. 1.e-2_8) then
            indform = 1
         elseif (plev(1,k) .lt. 1.e-1_8) then
            indform = 2
         elseif (plev(1,k) .lt. 1._8) then
            indform = 3
         elseif (plev(1,k) .lt. 10._8) then
            indform = 4
         elseif (plev(1,k) .lt. 100._8) then
            indform = 5
         elseif (plev(1,k) .lt. 1000._8) then
            indform = 6
         else
            indform = 7
         endif
         write(*,outform(indform)) k-1, plev(1,k), uflx(1,k), &
              difdflux(1,k), dirdflux(1,k), dflx(1,k), dflx(1,k)-uflx(1,k), hr(1,k)  ! Output are not the same fields...
      enddo
      write(*,9903)page
               
 9900 format(1x,'LEVEL PRESSURE   UPWARD FLUX   DIFDOWN FLUX  DIRDOWN FL&
     &UX  DOWNWARD FLUX   NET FLUX    HEATING RATE')
 9901 format(1x,'         mb          W/m2          W/m2          W/m2&
     &        W/m2          W/m2       degree/day')
 9902 format(1x,i3,3x,f11.6,4x,1p,2(g12.6,2x),g13.6,3x,g16.9,0p)
 9903 format(a)     
      ! ___ COMPLETED ABOVE ___ !



end program rrtmg_sw_test