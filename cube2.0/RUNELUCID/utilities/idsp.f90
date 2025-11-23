program inverse_displacement_field
  use parameters
  implicit none
  ! integer,parameter :: p=2 ! 1 cic ,2 is TSC order
  ! integer,parameter :: ndim=3
  integer(8) i,j,k,l,i_dim,iq(3),pid8,itx,ity,itz,nlast,ip,np
  integer cur_checkpoint,idx(ndim,p+1),k0,j0,i0
  integer(4),allocatable :: rhoc(:,:,:,:,:,:)
  real qpos(3),xpos(3),spos(3),dpos(3),zshift,dx(ndim,p+1),l3(ndim)
  real,allocatable :: vc(:,:,:,:,:,:,:),xtoq(:,:,:,:),stoq(:,:,:,:),sltoq(:,:,:,:),rhox(:,:,:),rhos(:,:,:),rhosl(:,:,:)
  integer(izipx),allocatable :: xp(:,:)
  integer(izipv),allocatable :: vp(:,:)
  integer(8),allocatable :: pid(:)

  call geometry
  print*, 'Inverse displacement field analysis on resolution:'
  print*, 'ng=',ng
  print*, 'checkpoint at:'
  open(11,file='../z_checkpoint.txt',status='old')
  
  do i=1,nmax_redshift
    read(11,end=71,fmt='(f8.4)') z_checkpoint(i)
    print*, z_checkpoint(i)
  enddo
  71 n_checkpoint=i-1
  close(11)
  allocate(rhoc(nt,nt,nt,nnt,nnt,nnt))
  allocate(xtoq(ndim,ng,ng,ng),stoq(ndim,ng,ng,ng),sltoq(ndim,ng,ng,ng),rhox(ng,ng,ng),rhos(ng,ng,ng),rhosl(ng,ng,ng))
  do cur_checkpoint= n_checkpoint-6,n_checkpoint-1
    sim%cur_checkpoint=cur_checkpoint
    print*, 'Start analyzing redshift ',z2str(z_checkpoint(cur_checkpoint))
    open(11,file=output_name('info'),status='old',action='read',access='stream')
    read(11) sim
    close(11)
    print*, 'nplocal =',sim%nplocal
    allocate(xp(3,sim%nplocal),vp(3,sim%nplocal),vc(3,nt,nt,nt,nnt,nnt,nnt),pid(sim%nplocal))
    open(11,file=output_name('xp'),status='old',action='read',access='stream')
    read(11) xp
    close(11)
    !open(11,file=output_name('vp'),status='old',action='read',access='stream')
    !read(11) vp
    !close(11)
    open(11,file=output_name('np'),status='old',action='read',access='stream')
    read(11) rhoc
    close(11)
    !open(11,file=output_name('vc'),status='old',action='read',access='stream')
    !read(11) vc
    !close(11)
    open(11,file=output_name('id'),status='old',action='read',access='stream')
    print*, output_name('id')
    read(11) pid
    close(11)
    print*,'check PID range:',minval(pid),maxval(pid)

    
    xtoq=0; stoq=0; sltoq=0; rhox=0; rhos=0; rhosl=0; nlast=0
    do itz=1,nnt
    do ity=1,nnt
    do itx=1,nnt
      do k=1,nt
      do j=1,nt
      do i=1,nt
        np=rhoc(i,j,k,itx,ity,itz)
        do l=1,np
          ip=nlast+l
          pid8=pid(ip)-1
          iq(3)=pid8/int(ng,4)**2     !this is tsc cic ?!
          iq(2)=(pid8-iq(3)*int(ng,4)**2)/int(ng,4)
          iq(1)=modulo(pid8,int(ng,4))

          qpos=iq+0.5 ! Lagrangian postion q
          xpos=(nt*([itx,ity,itz]-1)+[i,j,k]-1+(int(xp(:,ip)+ishift,izipx)+rshift)*x_resolution) * real(ng)/real(nc) ! Eulerian x
          !spos=xpos
          !zshift=vc(zdim,i,j,k,itx,ity,itz)+tan((pi*real(vp(zdim,ip)))/real(nvbin-1))/(sqrt(pi/2)/(sim%sigma_vi*vrel_boost))
          !zshift=zshift*sim%vsim2phys/sim%a/100 ! convert to km/sec and divide by 100a, in Mpc/h 
          !zshift=zshift*(ng/box) ! convert to find grid
          !zshift=zshift*1! 调节红移畸变系数 0.5
          !spos(zdim)=modulo(spos(zdim)+zshift,real(ng)); ! redshift space position s

          !slpos=xpos
          !z_lshift=vc(zdim,i,j,k,itx,ity,itz)+tan((pi*real(vp(zdim,ip)))/real(nvbin-1))/(sqrt(pi/2)/(sim%sigma_vi*vrel_boost))
          !zshift=zshift*sim%vsim2phys/sim%a/100
          !zshift=zshift*(ng/box) ! convert to find grid
          !zshift=zshift*1
          !slpos(zdim)=modulo(xpos(zdim)+z_lshift,real(ng)); !! redshift space linear position sl
          !open(15,file=output_name('v_L'),access='stream')
          !read(15) vl
          !close(15)

          ! x to q
          dpos=modulo(qpos-xpos+ng/2,real(ng))-ng/2
          idx(:,2)=floor(xpos)+1; idx(:,1)=idx(:,2)-1; idx(:,3)=idx(:,2)+1; idx=modulo(idx-1,ng)+1
          l3=xpos-floor(xpos); dx(:,3)=l3**2/2; dx(:,1)=(1-l3)**2/2; dx(:,2)=1-dx(:,1)-dx(:,3)
          do k0=1,p+1
          do j0=1,p+1
          do i0=1,p+1
            rhox(idx(1,i0),idx(2,j0),idx(3,k0))=rhox(idx(1,i0),idx(2,j0),idx(3,k0))+dx(1,i0)*dx(2,j0)*dx(3,k0)
            xtoq(:,idx(1,i0),idx(2,j0),idx(3,k0))=xtoq(:,idx(1,i0),idx(2,j0),idx(3,k0))+dx(1,i0)*dx(2,j0)*dx(3,k0)*dpos
          enddo
          enddo
          enddo

          ! s to q
          !dpos=modulo(qpos-spos+ng/2,real(ng))-ng/2
          !idx(:,2)=floor(spos)+1; idx(:,1)=idx(:,2)-1; idx(:,3)=idx(:,2)+1; idx=modulo(idx-1,ng)+1
          !l3=spos-floor(spos); dx(:,3)=l3**2/2; dx(:,1)=(1-l3)**2/2; dx(:,2)=1-dx(:,1)-dx(:,3)
          !do k0=1,p+1
          !do j0=1,p+1
          !do i0=1,p+1
          !  rhos(idx(1,i0),idx(2,j0),idx(3,k0))=rhos(idx(1,i0),idx(2,j0),idx(3,k0))+dx(1,i0)*dx(2,j0)*dx(3,k0)
          !  stoq(:,idx(1,i0),idx(2,j0),idx(3,k0))=stoq(:,idx(1,i0),idx(2,j0),idx(3,k0))+dx(1,i0)*dx(2,j0)*dx(3,k0)*dpos
          !enddo
          !enddo
          !enddo
          ! slinear to q
          !dpos=modulo(qpos-slpos+ng/2,real(ng))-ng/2
          !idx(:,2)=floor(spos)+1; idx(:,1)=idx(:,2)-1; idx(:,3)=idx(:,2)+1; idx=modulo(idx-1,ng)+1 !interpletion
          !l3=slpos-floor(slpos); dx(:,3)=l3**2/2; dx(:,1)=(1-l3)**2/2; dx(:,2)=1-dx(:,1)-dx(:,3)
          !do k0=1,p+1
          !do j0=1,p+1
          !do i0=1,p+1
            !rhosl(idx(1,i0),idx(2,j0),idx(3,k0))=rhosl(idx(1,i0),idx(2,j0),idx(3,k0))+dx(1,i0)*dx(2,j0)*dx(3,k0)
            !sltoq(:,idx(1,i0),idx(2,j0),idx(3,k0))=sltoq(:,idx(1,i0),idx(2,j0),idx(3,k0))+dx(1,i0)*dx(2,j0)*dx(3,k0)*dpos
          !enddo
          !enddo
          !enddo
        enddo
        nlast=nlast+np
      enddo
      enddo
      enddo
    enddo
    enddo
    enddo
    deallocate(xp,vp,vc,pid)

    do i_dim=1,3
      xtoq(i_dim,:,:,:)=xtoq(i_dim,:,:,:)/(rhox+0.0001)
      !stoq(i_dim,:,:,:)=stoq(i_dim,:,:,:)/(rhos+0.0001)
      print*, '  dimension',int(i_dim,1),'min,max values ='
      print*, '  xtoq',minval(xtoq(i_dim,:,:,:)), maxval(xtoq(i_dim,:,:,:))
      !print*, '  stoq',minval(stoq(i_dim,:,:,:)), maxval(stoq(i_dim,:,:,:))
    enddo

    print*,'Write idsp into file'
    print*,output_name('xtoq')
    open(11,file=output_name('xtoq'),access='stream')
    write(11) xtoq
    close(11)
    !open(11,file=output_name('stoq'),access='stream')
    !write(11) stoq
    !close(11)
    
  enddo ! do cur_checkpoint= n_checkpoint,n_checkpoint
  deallocate(rhox,rhos,xtoq,stoq)
  deallocate(rhoc)

  print*,'idisplacement done'
  
end
