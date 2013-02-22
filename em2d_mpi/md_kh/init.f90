module init

  use const
  use mpi_set
  use random_gen

  implicit none

  private

  public :: init__set_param

  integer, public, parameter   :: nroot=0
  integer, allocatable, public :: np2(:,:)
  integer, public              :: itmax, it0, intvl1, intvl2
  real(8), public              :: delx, delt, gfac
  real(8), public              :: c, q(nsp), r(nsp)
  real(8), allocatable, public :: uf(:,:,:)
  real(8), allocatable, public :: up(:,:,:,:)
  real(8), allocatable, public :: gp(:,:,:,:)
  character(len=64),  public    :: dir
  character(len=64), public    :: file12
  real(8)                      :: pi, n0, v0, b0, x0, vti, vte, rtemp, delv, dr, br, theta_b


contains

  
  subroutine init__set_param

    use fio, only : fio__input, fio__param
    real(8) :: fgi, fpi, alpha, beta, va, fpe, fge, rgi, rge, ldb
    character(len=64) :: file9 
    character(len=64) :: file11

!************** MPI settings  *******************!
    call mpi_set__init(nxgs,nxge,nygs,nyge,nproc)

    allocate(np2(nys:nye,nsp))
    allocate(uf(6,nxs-2:nxe+2,nys-2:nye+2))
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
    itmax  = 27000
    intvl1 = 900
    intvl2 = 900
    dir    = './dat/'
    file9  = 'init_param.dat'
    file12 = 'energy.dat'
    gfac   = 0.501

    it0    = 0
    if(it0 /= 0)then
       !start from the past calculation
       write(file11,'(i6.6,a,i3.3,a)') it0,'_rank=',nrank,'.dat'
       call fio__input(up,uf,np2,c,q,r,delt,delx,it0,                             &
                       np,nxgs,nxge,nygs,nyge,nxs,nxe,nys,nye,nsp,bc,nproc,nrank, &
                       dir,file11)
       return
    endif

!*********************************************************************
!   r(1)  : ion mass             r(2)  : electron mass
!   q(1)  : ion charge           q(2)  : electron charge
!   c     : speed of light       ldb   : debye length
!
!   rgi   : ion Larmor radius    rge   : electron Larmor radius
!   fgi   : ion gyro-frequency   fge   : electron gyro-frequency
!   vti   : ion thermal speed    vte   : electron thermal speed
!   b0    : magnetic field       
!  
!   alpha : wpe/wge
!   beta  : ion plasma beta
!   rtemp : Te/Ti
!*********************************************************************
    delx = 1.0
    c    = 1.0
    delt = 1.0*delx/c
    ldb  = delx*1.0

    r(1) = 16.0
    r(2) = 1.0

    alpha = 2.0
    beta  = 1.0
    rtemp = 1.0

    fpe = dsqrt(beta*rtemp)*c/(dsqrt(2.D0)*alpha*ldb)
    fge = fpe/alpha

    va  = fge/fpe*c*dsqrt(r(2)/r(1))
    rge = alpha*ldb*dsqrt(2.D0)
    rgi = rge*dsqrt(r(1)/r(2))/dsqrt(rtemp)
    vte = rge*fge
    vti = vte*dsqrt(r(2)/r(1))/dsqrt(rtemp)
    v0  = 1.0*dsqrt(1.+5./6.*beta*(1.+rtemp))*va

    fgi = fge*r(2)/r(1)
    fpi = fpe*dsqrt(r(2)/r(1))

    !average number density at x=nxge
    n0 = 100

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
    q(1) = fpi*dsqrt(r(1)/(4.0*pi*n0))
    q(2) = -q(1)

    !Magnetic field strength
    b0 = fgi*r(1)*c/q(1)
    !Magnetic field strength ratio
    br = 1.0
    !Magnetic field orientation
    theta_b = 0.D0/180.D0*pi

    !position of the shear layer
    x0 = 0.5*(nxge+nxgs)*delx 
    !half shear width
    delv = rgi*2.
    
    call random_gen__init(nrank)
    call init__loading
    call fio__param(np,nsp,np2,                             &
                    nxgs,nxge,nygs,nyge,nys,nye,            &
                    c,q,r,n0,0.5*r(1)*vti**2,rtemp,fpe,fge, &
                    ldb,delt,delx,dir,file9,                &
                    nroot,nrank)

  end subroutine init__set_param


  subroutine init__loading

    use boundary, only : boundary__particle, boundary__field

    integer, parameter   :: nbin=nx*20, nchart=nbin, nbin_min=nbin/20
    integer              :: i, j, k, l, ii, isp, isn
    real(8), parameter   :: eps=1d-2
    real(8)              :: dn, dnmax
    real(8)              :: sd, r1, r2, left, right, dx, intf0
    real(8)              :: e0, tp, vfnc, efnc, bfnc, dfnc, gamp, tmp, vd
    real(8)              :: intf(nchart), cx(nchart)
    real(8)              :: bfabs, cdf_n, dcn, vy, gam0, dxvy, dvx, dvy, den, x, y

    !*** functions ***!
    !magnetic field strength
    bfabs(x) = +0.5*b0*((1.+br)+(1.-br)*tanh((x-x0)/delv))

    !velocity profile
    vy(x) = +0.5*v0*tanh((x-x0)/delv)

    !Lorentz factor of mean flow velocity
    gam0(x) = 1./dsqrt(1.-(vy(x)/c)**2)

    !gradient of velocity profile
    dxvy(x) = +0.5*v0/delv*cosh((x-x0)/delv)**(-2)

    !perturbation in vx
    dvx(x,y) = -0.05*v0*sin(2.*pi*(y-nygs)/(nyge-nygs+1))/cosh((x-x0)/delv)**2

    !number density profile
    den(x) = 0.5*n0*((1.+dr)+(1.-dr)*tanh((x-x0)/delv))

    !cdf of density profile
    cdf_n(x) = (+((1.+dr)*(   x)+(1.-dr)*delv*log(cosh((   x-x0)/delv))) &
                -((1.+dr)*(nxgs)+(1.-dr)*delv*log(cosh((nxgs-x0)/delv))) &
               )/((1.+dr)*(nxge-nxgs))  

    !charge density profile
    dcn(x) = abs(bfabs(x)*cos(theta_b)*gam0(x)**2/(4.*pi*q(1)*c)*dxvy(x))

    !drift velocity of ion
    dvy(x) = -c/(4.*pi*q(1)*den(x)*(rtemp+1.)) &
            *0.5*b0*cos(theta_b)/delv*(1.-br)/cosh((x-x0)/delv)**2

    !*** end of ***!

    !*** setting of fields ***!
    !magnetic field
    do j=nys,nye
    do i=nxs,nxe+bc
       uf(1,i,j) = 0.0D0
    enddo
    enddo
    do j=nys,nye
    do i=nxs,nxe
       uf(2,i,j) = bfabs(i*delx)*sin(theta_b)
       uf(3,i,j) = bfabs(i*delx)*cos(theta_b)
    enddo
    enddo

    !electric field
    do j=nys,nye
    do i=nxs,nxe
       uf(4,i,j) = -vy(i*delx)*uf(3,i,j)
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
    !*** end of ***!

    !additional charge to satisfy div(-Vy x Bz) = 4\pi(ni-ne)
    dn = (-uf(4,nxgs,nys)+uf(4,nxge,nys))/(4.*pi*q(1))
    dnmax = 0.0
    do i=nxgs,nxge
       dnmax = max(dnmax,abs(dcn(i*delx)))
    enddo
    dnmax = 0.5*dnmax

    !particle position - background particles
    do j=nys,nye
       do ii=1,np2(j,1)-dnmax*(nxge+bc-nxgs+1)*delx
          call random_number(r1)
          left  = nxgs*delx
          right = nxge*delx
          do while(right-left > eps)
             dx = (right+left)/2.
             if(cdf_n(dx) >= r1)then
                right = dx
             else
                left = dx
             endif
          enddo
          !x position
          up(1,ii,j,1) = dx
          up(1,ii,j,2) = up(1,ii,j,1)
          
          !y position
          call random_number(r1)
          up(2,ii,j,1) = dble(j)*delx+r1*delx
          up(2,ii,j,2) = up(2,ii,j,1)
       enddo
    enddo

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
             call random_number(r1)
             left  = 1
             right = nchart
             do while(right-left > 2.)
                if(intf((left+right)/2.+0.5) > r1)then
                   right = (left+right)/2.
                else
                   left = (left+right)/2.
                endif
             enddo
             up(1,ii,j,isp) = cx((left+right)/2.+0.5)

             call random_number(r1)
             up(2,ii,j,isp) = dble(j)*delx+delx*r1
          enddo
       enddo
       np2(nys:nye,isp) = np2(nys:nye,isp)-isn*dn/2.
    enddo

    !velocity
    !total pressure at x=nxge
    e0 = -0.5*v0*b0/c
    tp = 0.5*n0*(r(1)*vti**2+r(2)*vte**2)+(b0*b0-e0*e0)/(8.*pi)

    do isp=1,nsp
       do j=nys,nye
          do ii=1,np2(j,isp)
             vfnc = vy(up(1,ii,j,isp))
             dfnc = den(up(1,ii,j,isp))
             bfnc = bfabs(up(1,ii,j,isp))
             efnc = -bfnc*vfnc/c
             
             sd = dsqrt( 2.*(tp-(bfnc**2-efnc**2)/(8.*pi)) &
                        /(dfnc*r(1)*(1.+rtemp)) )/dsqrt(2.D0)

             if(isp == 2) then
                sd = sd*dsqrt(r(1)/r(2)*rtemp)
             endif

             !non-relativistic Maxwellian
             call random_gen__bm(r1,r2)
             up(3,ii,j,isp) = sd*r1
             up(4,ii,j,isp) = sd*r2

             call random_gen__bm(r1,r2)
             up(5,ii,j,isp) = sd*r1

             gamp = dsqrt(1.D0+(up(3,ii,j,isp)**2+up(4,ii,j,isp)**2+up(5,ii,j,isp)**2)/(c*c))
             vd = dvy(up(1,ii,j,isp))
             if(isp == 2) vd = -vd*rtemp

             call random_number(r1)
             if(up(3,ii,j,isp)*dvx(up(1,ii,j,isp),up(2,ii,j,isp)) >= 0.)then
                up(3,ii,j,isp) = up(3,ii,j,isp)+dvx(up(1,ii,j,isp),up(2,ii,j,isp))*gamp
             else
                if(r1 < (-dvx(up(1,ii,j,isp),up(2,ii,j,isp))*up(3,ii,j,isp)/gamp))then
                   up(3,ii,j,isp) = -up(3,ii,j,isp)+dvx(up(1,ii,j,isp),up(2,ii,j,isp))*gamp
                else
                   up(3,ii,j,isp) = up(3,ii,j,isp)+dvx(up(1,ii,j,isp),up(2,ii,j,isp))*gamp
                endif
             endif
             call random_number(r1)
             if(up(4,ii,j,isp)*(+vy(up(1,ii,j,isp))+vd) >= 0.)then
                up(4,ii,j,isp) = +up(4,ii,j,isp)+(+vy(up(1,ii,j,isp))+vd)*gamp
             else
                if(r1 < (-(+vy(up(1,ii,j,isp))+vd)*up(4,ii,j,isp)/gamp))then
                   up(4,ii,j,isp) = -up(4,ii,j,isp)+(+vy(up(1,ii,j,isp))+vd)*gamp
                else
                   up(4,ii,j,isp) = +up(4,ii,j,isp)+(+vy(up(1,ii,j,isp))+vd)*gamp
                endif
             endif

          enddo
       enddo
    enddo

    call boundary__particle(up,                                        &
                            np,nsp,np2,nxgs,nxge,nygs,nyge,nys,nye,bc, &
                            nup,ndown,nstat,mnpi,mnpr,ncomw,nerr)

  end subroutine init__loading


end module init
