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
  character(len=64), public    :: dir
  character(len=64), public    :: file12
  real(8)                      :: pi, v0, x0, y0, vti, vte, rtemp, dr, br, theta_b
  real(8) :: ncs, nbg  ! Harris sheet density and the background density
  real(8) :: vdi, vde  ! drift speed
  real(8) :: b0   ! B0 for harris fields
  real(8) :: lcs  ! current sheet half thickness


contains

  
  subroutine init__set_param

    use fio, only : fio__input, fio__param
    real(8) :: fgi, fpi, alpha, va, fpe, fge, rgi, rge
    real(8) :: ldb       ! Debye length
    character(len=64) :: file9 
    character(len=64) :: file11

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
    pi     = 4.d0*atan(1.0)
    itmax  = 20000
    intvl1 = 500
    intvl2 = 50
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
!   c     : speed of light
!   ldb   : debye length         lcs   : current sheet width
!
!   rgi   : ion Larmor radius    rge   : electron Larmor radius
!   fgi   : ion gyro-frequency   fge   : electron gyro-frequency
!   vti   : ion thermal speed    vte   : electron thermal speed
!   vdi   : ion drift speed      vde   : electron drift speed
!   b0    : magnetic field       
!  
!   alpha : wpe/wge
!   rtemp : Te/Ti
!*********************************************************************
!    Reconnection setup  (by S. Zenitani, NAOJ)   2012/7/19
!*********************************************************************
    delx = 1.0
    c    = 1.0
    delt = 0.5*delx/c
    ldb  = delx

    r(1) = 16.d0
    r(2) =  1.d0

    alpha = 2.d0
    rtemp = 0.25d0

    vte = dsqrt(rtemp)*c/(dsqrt(1+rtemp)*alpha)
    vti = vte*dsqrt(r(2)/r(1))/dsqrt(rtemp)
    fpe = vte/ldb/dsqrt(2.d0)
    fpi = fpe*dsqrt(r(2)/r(1))
    fge = fpe/alpha
    fgi = fge*r(2)/r(1)
    va  = fge/fpe*c*dsqrt(r(2)/r(1))
    rge = vte/fge
    rgi = vti/fgi

    ! current sheet half thickness (fixed to 0.5 d_i)
    lcs = 0.5d0 * c/fpi
    ! current sheet density
    ncs = 250
    ! background density
    nbg = 50

    ! charge
    q(1) = fpi*dsqrt(r(1)/(4.d0*pi*ncs))
    q(2) = -q(1)  ! -fpe*dsqrt(r(2)/(4.d0*pi*ncs))
    ! Magnetic field amplitude
    b0  = fgi*r(1)*c/q(1)
    vdi = 1.d0 / ( 1.d0+rtemp ) * b0/(4*pi*lcs*q(1)*ncs)
    vde = -rtemp * vdi ! rtemp / (1.d0 + rtemp ) * b0 / ( 8*pi*lcs*q(1)*ncs )
    
    ! position of the X-point
    x0 = 0.5*(nxge+nxgs)*delx 
    y0 = 0.5*(nyge-nygs)*delx
    ! number of particles in each cell in y
    np2(nys:nye,1:nsp) = ceiling( nbg*(nxge-nxgs)*delx + ncs*2*lcs + 1d-6)

    if(nrank == nroot)then
       if(np2(nys,1) > np)then
          write(*,*)'Too large number of particles'
          stop
       endif
    endif
    
    call random_gen__init(nrank)
    call init__loading
    call fio__param(np,nsp,np2,                             &
                    nxgs,nxge,nygs,nyge,nys,nye,            &
                    c,q,r,ncs,0.5*r(1)*vti**2,rtemp,fpe,fge, &
                    ldb,delt,delx,dir,file9,                &
                    nroot,nrank)

  end subroutine init__set_param


  subroutine init__loading

    use boundary, only : boundary__particle, boundary__field

    integer              :: i, j, ii, ibg
    real(8)              :: r1, r2
    real(8)              :: b_harris, bx_pert, by_pert, jz, density
    real(8)              :: x, y
    real(8)              :: f1, f2
    real(8), parameter   :: e1 = 0.12d0  ! B1/B0: I recommend ~< O(0.5*nbg/ncs).

    ! ---------------- Utility functions ------------------
    ! magnetic field strength
    b_harris(x) = b0 * tanh((x-x0)/lcs)
    ! localized perturbation
    bx_pert(x,y) = +e1*b0 *((y-y0)/lcs) * exp(-((x-x0)**2+(y-y0)**2)/(2*lcs)**2)
    by_pert(x,y) = -e1*b0 *((x-x0)/lcs) * exp(-((x-x0)**2+(y-y0)**2)/(2*lcs)**2)
    ! density
    density(x) = ncs * cosh((x-x0)/lcs)**(-2) + nbg
    ! current jz_0 > 0, while jz_1 < 0
    jz(x,y)    = &
         +     b0/(4*pi*lcs) * cosh((x-x0)/lcs)**(-2) &
         -2*e1*b0/(4*pi*lcs) * ( 1.d0-((x-x0)**2+(y-y0)**2)/(2*lcs)**2 ) * exp(-((x-x0)**2+(y-y0)**2)/(2*lcs)**2)

    ! ---------------- Electromagnetic field ------------------
    do j=nys,nye
    do i=nxs,nxe+bc
       uf(1,i,j) = 0.0D0
    enddo
    enddo
    do j=nys,nye
    do i=nxs,nxe
       uf(2,i,j) = b_harris(i*delx)
       uf(3,i,j) = 0.d0
       uf(4,i,j) = 0.d0
    enddo
    enddo
    do j=nys,nye
    do i=nxs,nxe+bc
       uf(1,i,j) = uf(1,i,j) + bx_pert(i*delx,j*delx)
    enddo
    enddo
    do j=nys,nye
    do i=nxs,nxe
       uf(2,i,j) = uf(2,i,j) + by_pert(i*delx,j*delx)
    enddo
    enddo
    do j=nys,nye
    do i=nxs,nxe+bc
       uf(5,i,j) = 0.d0
       uf(6,i,j) = 0.d0
    enddo
    enddo

    call boundary__field(uf,                 &
                         nxs,nxe,nys,nye,bc, &
                         nup,ndown,mnpr,nstat,ncomw,nerr)
    ! ---------------- Electromagnetic field ------------------

    ! ---------------- Particles ------------------
    f1 =  1.d0 / ( (1.d0+rtemp) * q(1) )
    f2 = rtemp / ( (1.d0+rtemp) * q(2) )

    do j=nys,nye

       ibg = floor( nbg*(nxe+bc-nxs+1) + 1d-6)
       ! Background density: Uniform distribution
       do ii=1,ibg
          call random_number(r1)
          up(1,ii,j,1) = nxs*delx + r1*delx*(nxe+bc-nxs+1.d0)
          up(1,ii,j,2) = up(1,ii,j,1)
          call random_number(r1)
          up(2,ii,j,1) = dble(j)*delx + r1*delx
          up(2,ii,j,2) = up(2,ii,j,1)
       enddo
       ! Current sheet: 2nd order logistic distribution in x
       do ii=ibg+1,np2(j,1)
          call random_number(r1)
          up(1,ii,j,1) = lcs * 0.5d0*(log(r1)-log(1.d0-r1)) + x0
          up(1,ii,j,2) = up(1,ii,j,1)
          call random_number(r1)
          up(2,ii,j,1) = dble(j)*delx + r1*delx
          up(2,ii,j,2) = up(2,ii,j,1)
       enddo

       ! nonrelativistic Maxwellian 
       do ii=1,np2(j,1)
          call random_gen__bm(r1,r2)
          up(3,ii,j,1) = vti*r1
          up(4,ii,j,1) = vti*r2
          call random_gen__bm(r1,r2)
          up(5,ii,j,1) = vti*r1 + f1*jz(up(1,ii,j,1),up(2,ii,j,1))/density(up(1,ii,j,1))
          up(3,ii,j,2) = vte*r2
          call random_gen__bm(r1,r2)
          up(4,ii,j,2) = vte*r1
          up(5,ii,j,2) = vte*r2 + f2*jz(up(1,ii,j,2),up(2,ii,j,2))/density(up(1,ii,j,2))
       enddo
    enddo

    call boundary__particle(up,                                        &
                            np,nsp,np2,nxgs,nxge,nygs,nyge,nys,nye,bc, &
                            nup,ndown,nstat,mnpi,mnpr,ncomw,nerr)

  end subroutine init__loading


end module init
