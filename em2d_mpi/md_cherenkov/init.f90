module init

  use const
  use mpi_set
  use random_gen

  implicit none

  private

  public :: init__set_param

  integer, public, parameter :: nroot=0
  integer, allocatable, public :: np2(:,:)
  integer, public :: itmax, it0, intvl1, intvl2, intvl3
  real(8), public :: delx, delt, gfac
  real(8), public :: c
  real(8), allocatable, public :: uf(:,:,:)
  real(8), allocatable, public :: up(:,:,:,:)
  real(8), public :: q(nsp), r(nsp)
  real(8), allocatable, public :: gp(:,:,:,:)
  character(len=64), public :: dir
  character(len=64), public :: file12
  real(8)                   :: pi, vti, vte, va, v0, rtemp, fpe, fge, ldb, b0, n0
  real(8) :: gam0


contains

  
  subroutine init__set_param
    use fio, only : fio__input, fio__param
    real(8) :: fpi, n0
    real(8) :: L, rm, qw, q0(nsp)
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
!   intvl3  : printing interval for wave analysis
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
    itmax  = 600
    intvl1 = 100
    intvl2 = 5
    intvl3 = 5
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
!   L     : size of simulation box in electron skin depth
!
!   r(1)  : ion mass             r(2)  : electron mass
!   q(1)  : ion charge           q(2)  : electron charge
!   c     : speed of light       ldb   : debye length
!
!   fpi   : ion plasma frequency fpe   : electron plasma frequency
!   vti   : ion thermal speed    vte   : electron thermal speed
!   b0    : magnetic field       
!*********************************************************************
  
    delx = 1.0
    c    = 1.0
    delt = 0.5d0 * delx/c
    ldb  = delx
    !-- electron-positron plasma --
    ! positron
    r(1) = 1.0
    ! electron
    r(2) = 1.0
    
    ! thermal velocity
    vte = 0.1 * c
    vti = 0.1 * c
    ! Lorentz factor
    gam0 = 100
    v0   = c*dsqrt(1.-1./gam0**2)
    ! avarage number density at x = nxgs (lab.frame)
    n0 = 10

    ! number of particles
    np2(nys:nye,1) = n0*(nxe+bc-nxs+1)
    np2(nys:nye,2) = np2(nys:nye,1)

    if(nrank == nroot)then
       if(max(np2(nys,1), np2(nye,1), np) > np)then
          write(*,*)'Too large number of particles'
          stop
       endif
    endif

    ! charge
    if(nrank == nroot) n0 = dble(np2(nys,1))/dble((nxge-nxgs+1))
    call MPI_BCAST(n0,1,mnpr,nroot,ncomw,nerr)
    
    q(1) = vte/(ldb*dsqrt(2.0d0))* dsqrt(r(1) * gam0 / (4. * pi * n0))
    q(2) = -q(1)

    call random_gen__init(nrank)
    call init__loading
    call init__set_field
    call fio__param(np,nsp,np2,                     &
                    nxgs,nxge,nygs,nyge,nys,nye,    &
                    c,q,r,n0,0.5*r(1)*vti**2,rtemp,dsqrt(4.*pi*n0/gam0*q(1)**2/r(1)),q(1)*b0/(gam0*r(1)*c),    &
                    ldb,delt,delx,dir,file9,        &
                    nroot,nrank)

  end subroutine init__set_param

  
  
  subroutine init__loading

    use boundary, only : boundary__particle

    integer :: j, ii, isp
    real(8) :: sd, r1, r2

    !particle position
    isp=1
    do j=nys,nye
       do ii=1,np2(j,isp)
          call random_number(r1)
          up(1,ii,j,1) = nxs*delx+r1*delx*(nxe+bc-nxs+1.)
          up(1,ii,j,2) = up(1,ii,j,1)
          call random_number(r1)
          up(2,ii,j,1) = dble(j)*delx+delx*r1
          up(2,ii,j,2) = up(2,ii,j,1)
       enddo
    enddo

    !velocity
    !Maxwellian distribution
    do isp=1,nsp
       if(isp .eq. 1) then 
          sd = vti/dsqrt(2.0D0)
       endif
       if(isp .eq. 2) then
          sd = vte/dsqrt(2.0D0)
       endif

       do j=nys,nye
          do ii=1,np2(j,isp)
             call random_gen__bm(r1,r2)
             up(3,ii,j,isp) = sd*r1 + v0/sqrt(1-v0**2/c**2) 

             call random_gen__bm(r1,r2)
             up(4,ii,j,isp) = sd*r1
             up(5,ii,j,isp) = sd*r2
          enddo
       enddo
    enddo

    call boundary__particle(up,                                        &
                            np,nsp,np2,nxgs,nxge,nygs,nyge,nys,nye,bc, &
                            nup,ndown,nstat,mnpi,mnpr,ncomw,nerr)

  end subroutine init__loading


  
  subroutine init__set_field

    use boundary, only : boundary__field

    integer :: i, j

    !magnetic field
    do j=nys,nye
    do i=nxs,nxe+bc
       uf(1,i,j) = 0.0
    enddo
    enddo
    do j=nys,nye
    do i=nxs,nxe
       uf(2,i,j) = 0.0
       uf(3,i,j) = 0.0
    enddo
    enddo

    !electric field
    do j=nys,nye
    do i=nxs,nxe
       uf(4,i,j) = 0.0
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

  end subroutine init__set_field


end module init
