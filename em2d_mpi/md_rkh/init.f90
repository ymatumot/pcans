module init

  use const
  use mpi_set

  implicit none

  private

  public :: init__set_param

  integer, public, parameter :: nroot=0
  integer, allocatable, public :: np2(:,:)
  integer, public :: itmax, it0, intvl1, intvl2
  real(8), public :: delx, delt, gfac
  real(8), public :: c,q(nsp), r(nsp)
  real(8), allocatable, public :: uf(:,:,:)
  real(8), allocatable, public :: up(:,:,:,:)
  real(8), allocatable, public :: gp(:,:,:,:)
  character(len=6),  public :: dir
  character(len=10), public :: file12
  real(8)                   :: pi, n0, u0, b0, temp, rtemp, x0, delv, dr, br, theta_b


contains

  
  subroutine init__set_param

    use fio, only : fio__input, fio__param
    real(8) :: alpha, ldb, fpe, fge
    character(len=14) :: file9 
    character(len=19) :: file11

!************** MPI settings  *******************!
    call mpi_set__init(nxgs,nxge,nygs,nyge,nproc)

    allocate(np2(nys:nye,nsp))
    allocate(uf(6,nxs1:nxe1,nys1:nye1))
    allocate(up(5,np,nys:nye,nsp))
    allocate(gp(5,np,nys:nye,nsp))
!*********** End of MPI settings  ***************!

!*********************************************************************
!   time0   : start time (if time0 < 0, initial data from input.f)
!   itmax   : number of iteration
!   it0     : base count
!   intvl1  : storage interval for particles & fields
!   intvl2  : printing interval for energy variation
!   dir     : directory name for data output
!   file??  : output file name for unit number ??
!           :  9 - initial parameters
!           : 10 - for saving all data
!           : 11 - for starting from saved data
!           : 12 - for saving energy history
!   gfac    : implicit factor
!             gfac < 0.5 : unstable
!             gfac = 0.5 : no implicit
!             gfac = 1.0 : full implicit
!*********************************************************************
    pi     = 4.0*atan(1.0)
    itmax  = 0
    intvl1 = 400
    intvl2 = 200
    dir    = './dat/'
    file9  = 'init_param.dat'
    file12 = 'energy.dat'
    gfac   = 0.505

    it0    = 0
    if(it0 /= 0)then
       !start from the past calculation
       write(file11,'(a,i3.3,a)')'000100_rank=',nrank,'.dat'
       call fio__input(up,uf,c,q,r,delt,delx,it0,                                     &
                       np,np2,nxgs,nxge,nygs,nyge,nxs,nxe,nys,nye,nsp,bc,nproc,nrank, &
                       dir,file11)
       return
    endif

!*********************************************************************
!   r(1)  : ion mass             r(2)  : electron mass
!   q(1)  : ion charge           q(2)  : electron charge
!   c     : speed of light       ldb   : debye length
!
!   rge   : electron Larmor radius
!   fge   : electron gyro-frequency
!   b0    : magnetic field       
!  
!   alpha : wpe/wge
!   rtemp : Te/Ti
!*********************************************************************
    delx = 1.0
    c    = 1.0
    delt = 0.2
    ldb  = delx*4.0

    r(1) = 1.0
    r(2) = 1.0

    alpha = 1.0/dsqrt(2.D0)
    temp  = 3.*1.0*r(1)*c**2 ! proton temperature (T=3*Txx)
    rtemp = 1.0

    !flow speed
    u0  = 0.8*c

    !average number density at x=nxge (magnetosheath) in lab. frame
    n0 = 50.*dsqrt(1.+(0.5*u0/c)**2)

    !number density ratio of Nsph/Nsh
    dr = 1.0

    if(nrank == nroot)then
       if(0.5*(1.+dr)*n0*(nxe+bc-nxs+1) > np)then
          write(*,*)'Too large number of particles'
          stop
       endif
    endif

    !number of particles in each cell in y
    np2(nys:nye,1:nsp) = 0.5*n0*(1.+dr)*(nxe+bc-nxs+1)*delx

    !charge
    q(1) = dsqrt(temp/(4.*pi*n0*ldb**2))
    q(2) = -q(1)

    !plasma- gyro frequencies
    fpe = dsqrt(4.*pi*n0*q(1)**2/r(1))
    fge = fpe/alpha

    !Magnetic field strength
    b0 = fge*r(2)*c/abs(q(2))
    !Magnetic field strength ratio
    br = 1.0
    !Magnetic field orientation
    theta_b = 0.D0/180.D0*pi

    !position of the shear layer
    x0 = 0.5*(nxge+nxgs)*delx 

    !half shear width
    delv = (temp+r(1)*c**2)/(q(1)*b0)*4.0

    call init__loading_rela
    call fio__param(np,nsp,np2,                  &
                    nxgs,nxge,nygs,nyge,nys,nye, &
                    c,q,r,n0,temp,rtemp,fpe,fge, &
                    ldb,delt,delx,dir,file9,     &
                    nroot,nrank)

  end subroutine init__set_param


  subroutine init__loading_rela

    use boundary, only : boundary__particle, boundary__field

    integer, parameter :: nbin=nx*20, nchart=nbin, nbin_min=nbin/20
    integer :: i, j, ii, isp, k, l, n, isn
    integer, allocatable :: seed(:)
    real(8) :: dn, dnmax, left, right
    real(8) :: aa, bb, cc, dd, tmp, ut
    real(8) :: dx, intf0, gamp, u0y, v0y, dvx
    real(8) :: intf(nchart), cx(nchart)
    real(8) :: f, gf, u, tt
    real(8) :: uy, dxuy, dux, gam0, un, babs, dcn, x, y

    !*** functions ***!
    ! relativistic maxwell distribution function
    f(u,tt) = u**2*dexp(-dsqrt(1.0+(u/c)**2)*r(isp)*c**2/tt)

    ! gamma function
    gf(u,tt) = u**2*dexp(-r(isp)*c*u/tt)

    !velocity profile
    uy(x) = +0.5*u0*tanh((x-x0)/delv)

    !gradient of velocity profile
    dxuy(x) = +0.5*u0/delv*cosh((x-x0)/delv)**(-2)

    !Lorentz factor of mean flow velocity
    gam0(x) = dsqrt(1.+(uy(x)/c)**2)

    !perturbation in ux
    dux(x,y) = -0.05*u0*sin(2.*pi*(y-nygs)/(nyge-nygs+1))/cosh((x-x0)/delv)**2

    !number density profile in lab. frame
    un(x) = +0.5*n0*((1.+dr)+(1.-dr)*tanh((x-x0)/delv)) &
            *dsqrt(1.+(uy(x)/c)**2)

    !magnetic fields in rest frame
    babs(x) = +0.5*b0/dsqrt(1.+(0.5*u0/c)**2)*((1.+br)+(1.-br)*tanh((x-x0)/delv))

    !charge density profile
    dcn(x) = abs(babs(x)*cos(theta_b)/(4.*pi*q(1)*c)*dxuy(x))
    !*** end of ***!

    !*** field settings ***!
    !magnetic field
    do j=nys,nye
    do i=nxs,nxe+bc
       uf(1,i,j) = 0.0
    enddo
    enddo
    do j=nys,nye
    do i=nxs,nxe
       uf(2,i,j) = babs(i*delx)*sin(theta_b)
       uf(3,i,j) = babs(i*delx)*cos(theta_b)*gam0(i*delx)
    enddo
    enddo

    !electric field
    do j=nys,nye
    do i=nxs,nxe
       uf(4,i,j) = -uy(i*delx)/gam0(i*delx)/c*uf(3,i,j)
    enddo
    enddo
    do j=nys,nye
    do i=nxs,nxe+bc
       uf(5,i,j) = 0.0
       uf(6,i,j) = 0.0
    enddo
    enddo

    call boundary__field(uf,                 &
                         nxs,nxe,nys,nye,bc, &
                         nup,ndown,mnpr,nstat,ncomw,nerr)

    !*** particle loading ***!
    call random_seed(size=n)
    allocate(seed(n))
    seed(1:n) = nrank
    call random_seed(put=seed)
    deallocate(seed)
!!$    call random_seed()

    !Generation of random numbers for a number density profile;
    !CDF is evaluated by Simpson's rule
    dx = (nxge+bc-nxgs+1)*delx/nbin
    intf0 = un(nxgs*delx)
    do k=1,nbin-2,2
       intf0 = intf0+4.*un(nxgs*delx+dx*k)+2.*un(nxgs*delx+dx*(k+1))
    enddo
    k=nbin-1
    intf0 = intf0+4.*un(nxgs*delx+dx*k)+un(nxgs*delx+dx*(k+1))
    intf0 = intf0*dx/3.D0

    do l=1,nchart
       cx(l) = nxgs*delx+(nxge+bc-nxgs+1)*delx*(l-1)/(nchart-1)
       dx = (cx(l)-nxgs*delx)/(nbin_min+(nbin-nbin_min)*(l-1)/dble(nchart-1))
       intf(l) = un(nxgs*delx)
       do k=1,(nbin_min+(nbin-nbin_min)*(l-1)/dble(nchart-1))-2,2
          intf(l) = intf(l)+4.*un(nxgs*delx+dx*k)+2.*un(nxgs*delx+dx*(k+1))
       enddo
       k=(nbin_min+(nbin-nbin_min)*(l-1)/dble(nchart-1))-1
       intf(l) = intf(l)+4.*un(nxgs*delx+dx*k)+un(nxgs*delx+dx*(k+1))
       intf(l) = intf(l)*dx/3.D0
    enddo
    intf(1:nchart) = intf(1:nchart)/intf0

    !additional charge to satisfy div(-Vy/c x Bz) = 4\pi(ni-ne)
    dn = (-uf(4,nxgs,nys)+uf(4,nxge,nys))/(4.*pi*q(1))
    dnmax = 0.0
    do i=nxgs,nxge
       dnmax = max(dnmax,abs(dcn(i*delx)))
    enddo
    dnmax = 0.5*dnmax

    !particle position - background particles
    do j=nys,nye
       do ii=1,np2(j,1)-dnmax*(nxge+bc-nxgs+1)*delx
          call random_number(aa)
          left  = 1
          right = nchart
          do while(right-left > 2.)
             if(intf((left+right)/2.+0.5) > aa)then
                right = (left+right)/2.
             else
                left = (left+right)/2.
             endif
          enddo

          !x position
          up(1,ii,j,1) = cx((left+right)/2.+0.5)
          up(1,ii,j,2) = up(1,ii,j,1)

          !y position
          call random_number(bb)
          up(2,ii,j,1) = j*delx+bb*delx
          up(2,ii,j,2) = up(2,ii,j,1)
       enddo
    enddo

    !particle position - additional particles for charge correction
    do isp=1,nsp
       isn = sign(1.0D0,dn)
       if(isp == 2) isn = -isn

       !Generation of random numbers for a charge density profile;
       !CDF is evaluated by Simpson's rule
       dx = (nxge+bc-nxgs+1)*delx/nbin
       intf0 = isn*0.5*dcn(nxgs*delx)+dnmax
       do k=1,nbin-2,2
          intf0 = intf0+4.*(isn*0.5*dcn(nxgs*delx+dx*k)+dnmax) &
                       +2.*(isn*0.5*dcn(nxgs*delx+dx*(k+1))+dnmax)
       enddo
       k=nbin-1
       intf0 = intf0+4.*(isn*0.5*dcn(nxgs*delx+dx*k)+dnmax) &
                    +isn*0.5*dcn(nxgs*delx+dx*(k+1))+dnmax
       intf0 = intf0*dx/3.D0

       do l=1,nchart
          cx(l) = nxgs*delx+(nxge+bc-nxgs+1)*delx*(l-1)/(nchart-1)
          dx = (cx(l)-nxgs*delx)/(nbin_min+(nbin-nbin_min)*(l-1)/dble(nchart-1))
          intf(l) = isn*0.5*dcn(nxgs*delx)+dnmax
          do k=1,(nbin_min+(nbin-nbin_min)*(l-1)/dble(nchart-1))-2,2
             intf(l) = intf(l)+4.*(isn*0.5*dcn(nxgs*delx+dx*k)+dnmax) &
                              +2.*(isn*0.5*dcn(nxgs*delx+dx*(k+1))+dnmax)
          enddo
          k=(nbin_min+(nbin-nbin_min)*(l-1)/dble(nchart-1))-1
          intf(l) = intf(l)+4.*(isn*0.5*dcn(nxgs*delx+dx*k)+dnmax) &
                           +isn*0.5*dcn(nxgs*delx+dx*(k+1))+dnmax
          intf(l) = intf(l)*dx/3.D0
       enddo
       intf(1:nchart) = intf(1:nchart)/intf0

       do j=nys,nye
          do ii=np2(j,isp)-dnmax*(nxge+bc-nxgs+1)*delx+1,np2(j,isp)-isn*dn/2
             call random_number(aa)
             left  = 1
             right = nchart
             do while(right-left > 2.)
                if(intf((left+right)/2.+0.5) > aa)then
                   right = (left+right)/2.
                else
                   left = (left+right)/2.
                endif
             enddo
             up(1,ii,j,isp) = cx((left+right)/2.+0.5)

             call random_number(bb)
             up(2,ii,j,isp) = dble(j)*delx+delx*bb
          enddo
       enddo
       np2(nys:nye,isp) = np2(nys:nye,isp)-isn*dn/2.
    enddo

    !velocity
    !relativistic Maxwellian distribution
    do isp=1,nsp
       do j=nys,nye
          do ii=1,np2(j,isp)

             if(isp .eq. 1) then 
                tmp = temp
             endif
             if(isp .eq. 2) then
                tmp = temp*rtemp
             endif

             !acceptance-rejection method based on Gamma function effective for T/mc^2 >> 1
             call random_number(aa)
             call random_number(bb)
             call random_number(cc)
             ut = tmp/(r(isp)*c)*(-dlog(aa)-dlog(bb)-dlog(cc))
             call random_number(dd)

             do while(dd > f(ut,tmp)/gf(ut,tmp))
                call random_number(aa)
                call random_number(bb)
                call random_number(cc)
                ut = tmp/(r(isp)*c)*(-dlog(aa)-dlog(bb)-dlog(cc))
                call random_number(dd)
             enddo

             call random_number(bb)
             call random_number(cc)
             up(3,ii,j,isp) = ut*(2.*bb-1.)
             up(4,ii,j,isp) = ut*2.*dsqrt(bb*(1.-bb))*cos(2.*pi*cc)
             up(5,ii,j,isp) = ut*2.*dsqrt(bb*(1.-bb))*sin(2.*pi*cc)

             call random_number(dd)
             gamp = dsqrt(1.D0+ut*ut/(c*c))
             v0y  = uy(up(1,ii,j,isp))/gam0(up(1,ii,j,isp))
             dvx  = dux(up(1,ii,j,isp),up(2,ii,j,isp))/gam0(up(1,ii,j,isp))

             if(up(3,ii,j,isp)*dvx >= 0.)then
                up(3,ii,j,isp) = gam0(up(1,ii,j,isp))*(+up(3,ii,j,isp)+dvx*gamp)
             else
                if(dd < (-dvx*up(3,ii,j,isp)/gamp))then
                   up(3,ii,j,isp) = gam0(up(1,ii,j,isp))*(-up(3,ii,j,isp)+dvx*gamp)
                else
                   up(3,ii,j,isp) = gam0(up(1,ii,j,isp))*(+up(3,ii,j,isp)+dvx*gamp)
                endif
             endif
             call random_number(dd)
             if(up(4,ii,j,isp)*v0y >= 0.)then
                up(4,ii,j,isp) = gam0(up(1,ii,j,isp))*(+up(4,ii,j,isp)+v0y*gamp)
             else
                if(dd < (-v0y*up(4,ii,j,isp)/gamp))then
                   up(4,ii,j,isp) = gam0(up(1,ii,j,isp))*(-up(4,ii,j,isp)+v0y*gamp)
                else
                   up(4,ii,j,isp) = gam0(up(1,ii,j,isp))*(+up(4,ii,j,isp)+v0y*gamp)
                endif
             endif

          enddo
       enddo
    enddo

    call boundary__particle(up,                                        &
                            np,nsp,np2,nxgs,nxge,nygs,nyge,nys,nye,bc, &
                            nup,ndown,nstat,mnpi,mnpr,ncomw,nerr)

  end subroutine init__loading_rela



end module init
