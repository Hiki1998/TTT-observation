program spincorr_sdss
  use omp_lib
  !use parameters
  !use pencil_fft
  use iso_fortran_env, only : int64
  implicit none

  character(*),parameter :: mangadir='/mnt/18T/cube/usr_out/mingjiesheng/mangainfo/'
  character(*),parameter :: simudir='/mnt/18T/cube/usr_out/mingjiesheng/ELUCID_ratio2_PP/image1/'

  real,parameter :: box=500  ! simulation scale /dim, in unit of Mpc/h
  integer,parameter :: ngrid=500 ! ELUCID grid number
  integer(8),parameter :: ncore=1
  real,parameter :: pi=4*atan(1.)

  integer(4) t1,t2,t_rate,nhalo,ihalo,hgrid(3),irand,inew,ig(3),ii,jj,i,j,k
  integer(8) plan_fft_fine,plan_ifft_fine,idx1(3,3)
  real kr,kx,ky,kz,pow,r_filter,rm,cen(3),dphi,qpos(3),g7(3),l3(3),dx(3,3),dx1(3),dx2(3)

  real rho_f(ngrid+2,ngrid,ngrid),phil(ngrid+2,ngrid,ngrid)
  complex crho_f(ngrid/2+1,ngrid,ngrid),phi_k(ngrid/2+1,ngrid,ngrid)
  real phis(ngrid+2,ngrid,ngrid) ! backup small scale phi
  integer i0,j0,k0,i1,j1,k1,i2,j2,k2,n_rsmall,n_ratio,itemp,l1,iband,nband
  real tide_s(3,3),tide_l(3,3),torque(3,3)
  real spin(3,ngrid,ngrid,ngrid),beta(ngrid,ngrid,ngrid)
  integer,allocatable :: idx(:,:,:)
  real,allocatable :: gdata(:,:)
  real,allocatable :: theta(:,:),corr(:,:,:),r_small(:),ratio_scale(:),qdata(:,:)
  real,allocatable :: idsp(:,:,:,:),wgt(:,:,:),jt(:,:,:,:,:),betaTT(:,:,:,:)
  equivalence(rho_f,crho_f,phil)


  !call omp_set_num_threads(ncore)
  cen=[30,370,370] ! earth location
  dphi=-39.*pi/180.
  nband=1
  rm=2.0 ! mass bin ratio
  n_rsmall=75
  n_ratio=1
  allocate(r_small(n_rsmall),ratio_scale(n_ratio),corr(n_rsmall,n_ratio,nband))
  do i1=1,n_rsmall
    r_small(i1)=0.2+0.2*(i1-1)
  enddo
  do i1=1,n_ratio
    ratio_scale(i1)=1.1+0.2*(i1-1)
  enddo

  print*, ''
  print*, 'program spincorr_manga on single node'
  print*, '  on',ncore,' cores'
  print*, '  Resolution ngrid =', ngrid
  print*, '  Box size', box
  print*, '-----------------------------------------'

  print*,''
  print*,'creating FFT plans'
  call system_clock(t1,t_rate)
  call create_cubefft_plan
  call system_clock(t2,t_rate)
  print*, '  elapsed time =',real(t2-t1)/t_rate,'secs';

  print*,''
  print*,'reading SDSS catalog'
  open(11,file=mangadir//'manga_cat_Elliptical_group_rsd.bin',status='old',access='stream')
  read(11) nhalo
  print*,'  nhalo =',nhalo
  allocate(gdata(3,nhalo))
  read(11) gdata
  close(11)

  print*,''
  print*,'rotate position and spin'
  do ihalo=1,nhalo
    g7 = gdata(:,ihalo)
    gdata(3,ihalo) = cen(3) + g7(1)*cos(dphi)      + g7(2)*sin(dphi);
    gdata(2,ihalo) = cen(2) + g7(1)*cos(dphi+pi/2) + g7(2)*sin(dphi+pi/2);
    gdata(1,ihalo) = cen(1) + g7(3);
    !gdata(4,ihalo) =          g7(4);
    !gdata(7,ihalo) =          g7(5)*cos(dphi)      + g7(6)*sin(dphi);
    !gdata(6,ihalo) =          g7(5)*cos(dphi+pi/2) + g7(6)*sin(dphi+pi/2);
    !gdata(5,ihalo) =          g7(7);
  enddo

  print*,''
  print*,'reading weight file'
  allocate(wgt(ngrid+2,ngrid,ngrid))
  open(11,file=mangadir//'Weight5hR5_0500.bin',status='old',access='stream')
  read(11) wgt
  close(11)

  print*,''
  print*,'reading remapping field'
  allocate(idsp(3,ngrid,ngrid,ngrid))
  open(11,file=simudir//'0.000_xtoq_1.bin',status='old',access='stream')
  read(11) idsp
  close(11)

  print*,''
  print*,'remapping to q'
  allocate(qdata(6,nhalo))
  inew=0
  qdata=0
  do ihalo=1,nhalo
    ig=ceiling(gdata(1:3,ihalo))
    qpos=gdata(1:3,ihalo)
    if (minval(ig)>=1 .and. maxval(ig)<=ngrid) then
    if (wgt(ig(1),ig(2),ig(3))>0.5) then
      !qpos=gdata(1:3,ihalo)+idsp(:,ig(1),ig(2),ig(3)) ! q space NGP
      idx1(:,2) = ig
      idx1(:,1)=idx1(:,2)-1 
      idx1(:,3)=idx1(:,2)+1 
      idx1=modulo(idx1-1,ngrid)+1
      l3=(gdata(1:3,ihalo))-floor(gdata(1:3,ihalo))
      dx(:,1)=(1-l3)**2/2 
      dx(:,3)=l3**2/2
      dx(:,2)=1-dx(:,1)-dx(:,3)
      do k=1,3
      do j=1,3
      do i=1,3 ! q space TSC
        qpos=qpos+idsp(:,idx1(1,i),idx1(2,j),idx1(3,k))*dx(1,i)*dx(2,j)*dx(3,k);
      enddo
      enddo
      enddo
      if (minval(qpos)>0. .and. maxval(qpos)<real(ngrid) .and. gdata(5,ihalo)>-10000) then
        inew=inew+1
        qdata(1:3,inew)=gdata(1:3,ihalo)
        qdata(4:6,inew)=qpos
        !qdata(8:10,inew)=gdata(5:7,ihalo)
      else
        print*, 'error index ',ihalo 
      endif
    endif
    endif
  enddo
  deallocate(wgt,idsp,gdata)

  print*, '  converted',inew,' galaxies, min & max in qx qy qz ='
  print*, '  ',minval(qdata(4,1:inew)),maxval(qdata(4,1:inew))
  print*, '  ',minval(qdata(5,1:inew)),maxval(qdata(5,1:inew))
  print*, '  ',minval(qdata(6,1:inew)),maxval(qdata(6,1:inew))

  print*,''
  print*,'  writing to file'
  open(10,file=simudir//'qdata.bin',status='replace',access='stream')
  write(10) inew
  write(10) qdata(:,1:inew)
  close(10)

  print*,''
  print*,'generating grid index'
  nhalo=inew
  allocate(idx(3,nhalo,nband),theta(nhalo,nband))
  do ihalo=1,nhalo
    idx(:,ihalo,1)=ceiling(qdata(4:6,ihalo)) ! q space
    idx(:,ihalo,2)=ceiling(qdata(1:3,ihalo)) ! s space
     ! below are randomized positions to test the error
    !irand=modulo(ihalo+nhalo/2,nhalo)+1
    !idx(1,ihalo,3)=ceiling(qdata(5,irand)+ngrid/2)
    !idx(2,ihalo,3)=ceiling(qdata(6,irand)+ngrid/2)
    !idx(3,ihalo,3)=ceiling(qdata(4,irand)+ngrid/2)
    !idx(:,ihalo,3)=modulo(idx(:,ihalo,3),ngrid)+1
  enddo
  if (maxval(idx)>ngrid .or. minval(idx)<1) then
    print*, 'error: idx out of range'
    print*, minval(idx), maxval(idx)
    print*,minloc(idx)
    stop
  endif

  print*,''
  print*,'reading ELUCID IC'
  open(11,file='../cxyz_251_500_500.bin',status='old',access='stream')
  read(11) rho_f
  close(11)

  print*,''
  print*,'convert to potential'
  do k1=1,ngrid
  do j1=1,ngrid
  do i1=1,ngrid/2+1
    kz=mod(k1+ngrid/2-1,ngrid)-ngrid/2
    ky=mod(j1+ngrid/2-1,ngrid)-ngrid/2
    kx=i1-1
    kz=2*sin(pi*kz/ngrid)
    ky=2*sin(pi*ky/ngrid)
    kx=2*sin(pi*kx/ngrid)

    kr=kx**2+ky**2+kz**2
    kr=max(kr,1.0/ngrid**2)
    !kr=2*pi*kr/ngrid ! without sinc function
    pow=-4*pi/kr
    crho_f(i1,j1,k1)=crho_f(i1,j1,k1)*pow
  enddo
  enddo
  enddo
  crho_f(1,1,1)=0
  phi_k=crho_f ! backup Fourier potential in phi_k

  print*,''
  print*,'reconstructing and correlating spin'
  allocate(jt(3,nhalo,n_rsmall,n_ratio,nband))
  allocate(betaTT(nhalo,n_rsmall,n_ratio,nband))
  do ii=1,n_rsmall
    print*, '  progress',ii,'/',n_rsmall
    call gaussian_fourier_filter(phi_k,r_small(ii))
    call sfftw_execute(plan_ifft_fine)
    phis=rho_f
    do jj=1,n_ratio
      call system_clock(t1,t_rate) ! tic
      call gaussian_fourier_filter(phi_k,r_small(ii)*ratio_scale(jj))
      call sfftw_execute(plan_ifft_fine) ! phil is equivalenced to rho_f
      call spinfield
      do iband=1,nband
      do ihalo=1,nhalo
        jt(:,ihalo,ii,jj,iband)=spin(:,idx(1,ihalo,iband),idx(2,ihalo,iband),idx(3,ihalo,iband))
        betaTT(ihalo,ii,jj,iband)=beta(idx(1,ihalo,iband),idx(2,ihalo,iband),idx(3,ihalo,iband))
        !theta(ihalo,iband)=ccc(qdata(8:10,ihalo),spin(:,idx(1,ihalo,iband),idx(2,ihalo,iband),idx(3,ihalo,iband)))
        !print*,ihalo,qdata(8,ihalo)
      enddo
      enddo
      !do iband=1,nband
      !  corr(ii,jj,iband)=sum(theta(:,iband))/nhalo
      !enddo
      call system_clock(t2,t_rate) ! toc
      !print*, r_small(ii), ratio_scale(jj), corr(ii,jj,:), real(t2-t1)/t_rate,'secs';

    enddo
  enddo

  open(11,file=simudir//'result.bin',status='replace',access='stream')
  write(11) n_rsmall,n_ratio,r_small(1:n_rsmall),ratio_scale(1:n_ratio)
  !write(11) corr
  write(11) jt
  write(11) betaTT
  close(11)

call destroy_cubefft_plan


contains

  real function ccc(vec_i,vec2)
    implicit none
    real vec_i(3),vec2(3)
    vec_i=vec_i/norm2(vec_i)
    vec2=vec2/norm2(vec2)
    ccc=sum(vec_i*vec2)
  endfunction

  subroutine spinfield
    implicit none
    save

    do k1=1,ngrid
      k2=modulo(k1,ngrid)+1
      k0=modulo(k1-2,ngrid)+1
    do j1=1,ngrid
      j2=modulo(j1,ngrid)+1
      j0=modulo(j1-2,ngrid)+1
    do i1=1,ngrid
      i2=modulo(i1,ngrid)+1
      i0=modulo(i1-2,ngrid)+1
!print*,phis(1,1,1);stop
      tide_s(1,1)=phis(i2,j1,k1)-2*phis(i1,j1,k1)+phis(i0,j1,k1)
      tide_s(2,2)=phis(i1,j2,k1)-2*phis(i1,j1,k1)+phis(i1,j0,k1)
      tide_s(3,3)=phis(i1,j1,k2)-2*phis(i1,j1,k1)+phis(i1,j1,k0)
      tide_s(1,2)=(phis(i2,j2,k1)+phis(i0,j0,k1)-phis(i2,j0,k1)-phis(i0,j2,k1))/4
      tide_s(2,3)=(phis(i1,j2,k2)+phis(i1,j0,k0)-phis(i1,j2,k0)-phis(i1,j0,k2))/4
      tide_s(3,1)=(phis(i2,j1,k2)+phis(i0,j1,k0)-phis(i2,j1,k0)-phis(i0,j1,k2))/4
      tide_s(2,1)=tide_s(1,2)
      tide_s(3,2)=tide_s(2,3)
      tide_s(1,3)=tide_s(3,1)
!print*,tide_s
      tide_l(1,1)=phil(i2,j1,k1)-2*phil(i1,j1,k1)+phil(i0,j1,k1)
      tide_l(2,2)=phil(i1,j2,k1)-2*phil(i1,j1,k1)+phil(i1,j0,k1)
      tide_l(3,3)=phil(i1,j1,k2)-2*phil(i1,j1,k1)+phil(i1,j1,k0)
      tide_l(1,2)=(phil(i2,j2,k1)+phil(i0,j0,k1)-phil(i2,j0,k1)-phil(i0,j2,k1))/4
      tide_l(2,3)=(phil(i1,j2,k2)+phil(i1,j0,k0)-phil(i1,j2,k0)-phil(i1,j0,k2))/4
      tide_l(3,1)=(phil(i2,j1,k2)+phil(i0,j1,k0)-phil(i2,j1,k0)-phil(i0,j1,k2))/4
      tide_l(2,1)=tide_l(1,2)
      tide_l(3,2)=tide_l(2,3)
      tide_l(1,3)=tide_l(3,1)
!print*,tide_l
!stop
      torque=-matmul(tide_s,tide_l)
      spin(1,i1,j1,k1)=-torque(2,3)+torque(3,2)
      spin(2,i1,j1,k1)=-torque(3,1)+torque(1,3)
      spin(3,i1,j1,k1)=-torque(1,2)+torque(2,1)
      beta(i1,j1,k1)=norm2(spin(:,i1,j1,k1))/sum(tide_s*tide_l)
      spin(:,i1,j1,k1)=spin(:,i1,j1,k1)/norm2(spin(:,i1,j1,k1))
    enddo
    enddo
    enddo

  endsubroutine

  subroutine gaussian_fourier_filter(phi_k,r_filter)
    ! apply Gaussian widxow function with radius r_filter to Fourier space phi_k
    ! returns Fourier space crho_f
    implicit none
    complex phi_k(ngrid/2+1,ngrid,ngrid)
    real r_filter
    integer i01,j01,k01
    crho_f=0
    do k01=1,ngrid
    do j01=1,ngrid
    do i01=1,ngrid/2+1
      kz=mod(k01+ngrid/2-1,ngrid)-ngrid/2
      ky=mod(j01+ngrid/2-1,ngrid)-ngrid/2
      kx=i01-1
      kr=sqrt(kx**2+ky**2+kz**2)
      kr=max(kr,1.0)
      kr=2*pi*kr/box
      pow=exp(-kr**2*r_filter**2/2)**0.25 ! apply E-mode widxow function
      crho_f(i01,j01,k01)=phi_k(i01,j01,k01)*pow
!if (k01==111.and.j01==111.and.i01==2) then
!  print*,kx,ky,kz,kr
!  print*,phi_k(i01,j01,k01)
!  print*,crho_f(i01,j01,k01)
!endif
    enddo
    enddo
    enddo
    crho_f(1,1,1)=0 ! DC frequency
  endsubroutine

  subroutine create_cubefft_plan
    use,intrinsic :: ISO_C_BINDING
    use omp_lib
    implicit none
    save
    !include 'fftw3.f'
    include 'fftw3.f03'
    integer istat,icore

#ifndef macbook
    call sfftw_init_threads(istat)
    print*, 'sfftw_init_threads status',istat
    icore=omp_get_max_threads()
    print*, 'omp_get_max_threads() =',icore
    !call sfftw_plan_with_nthreads(icore)
    call sfftw_plan_with_nthreads(64)
#endif

    call sfftw_plan_dft_r2c_3d(plan_fft_fine,ngrid,ngrid,ngrid,rho_f,rho_f,FFTW_MEASURE)
    call sfftw_plan_dft_c2r_3d(plan_ifft_fine,ngrid,ngrid,ngrid,rho_f,rho_f,FFTW_MEASURE)
  endsubroutine create_cubefft_plan

  subroutine destroy_cubefft_plan
    use,intrinsic :: ISO_C_BINDING
    implicit none
    save
    !include 'fftw3.f'
    include 'fftw3.f03'
    call sfftw_destroy_plan(plan_fft_fine)
    call sfftw_destroy_plan(plan_ifft_fine)
#ifndef macbook
    call fftw_cleanup_threads()
#endif
  endsubroutine destroy_cubefft_plan

endprogram
