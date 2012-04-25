module init

  use const
  use mpi_set
  use random_gen

  implicit none

  private

  public :: init__set_param

  integer, public, parameter :: nroot=0
  integer, allocatable, public :: np2(:,:)
  integer, public :: itmax, it0, intvl1, intvl2
  integer, public :: bcp        
  real(8), public :: delx, delt, gfac
  real(8), public :: c
  real(8), allocatable, public :: uf(:,:)
  real(8), allocatable, public :: up(:,:,:,:)
  real(8), public :: q(nsp), r(nsp)
  !gx, gv, are temporal spaces used for the time integration
  real(8), allocatable, public :: gp(:,:,:,:) !just for initialization
  character(len=64),  public :: dir
  character(len=64), public :: file12
  real(8)                   :: pi, vti, vte, va, rtemp, t_ani, fpe, fge, rgi, rge, ldb, b0


contains

  
  subroutine init__set_param

    use fio, only : fio__input, fio__param
    integer :: n0
    real(8) :: fgi, fpi, alpha, beta
    character(len=64) :: file9 
    character(len=64) :: file11

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
    itmax  = 50000
    intvl1 = 50
    intvl2 = 500
    dir    = './dat/'
    file9  = 'init_param.dat'
    file12 = 'energy.dat'
    gfac   = 0.505

!************** MPI settings  *******************!
    call mpi_set__init(nxgs,nxge,bc,nproc)
    if(nrank == nproc-1)then
       if(bc == -1) bcp = -1
       if(bc ==  0) bcp = 0
    else
       bcp = 0
    endif
    allocate(np2(nxs:nxe+bcp,nsp))
    allocate(uf(6,nxs1:nxe1))
    allocate(up(4,np,nxs:nxe+bcp,nsp))
    allocate(gp(4,np,nxs:nxe+bcp,nsp))
!*********** End of MPI settings  ***************!

    it0    = 0
    if(it0 /= 0)then
       !start from the past calculation
       write(file11,'(a,i3.3,a)')'002400_rank=',nrank,'.dat'
       call fio__input(up,uf,np2,c,q,r,delt,delx,it0, &
                       np,nxgs,nxge,nxs,nxe,nxs1,nxe1,nsp,bc,bcp,nproc,nrank,dir,file11)
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
!   t_ani : Te_perp/ Te_para
!*********************************************************************
    delx = 1.0
    c    = 1.0
    delt = 0.5
    ldb  = delx

    r(1) = 1837.0
    r(2) = 1.0

    alpha = 5.0
    beta  = 2.0
    rtemp = 1.0
    t_ani = 3.0

    fpe = dsqrt(beta*rtemp)*c/(dsqrt(2.D0)*alpha*ldb)
    fge = fpe/alpha

    va  = fge/fpe*c*dsqrt(r(2)/r(1))
    rge = fpe/fge*ldb*dsqrt(2.D0)
    rgi = rge*dsqrt(r(1)/r(2))/dsqrt(rtemp)

    vte = rge*fge
    vti = vte*dsqrt(r(2)/r(1))/dsqrt(rtemp)

    fgi = fge*r(2)/r(1)
    fpi = fpe*dsqrt(r(2)/r(1))

    np2(nxs:nxe+bcp,1) = 1000
    np2(nxs:nxe+bcp,2) = np2(nxs:nxe+bcp,1)

    if(nrank == nroot)then
       if(max(np2(nxs,1), np2(nxe+bcp,1), np) > np)then
          write(*,*)'Too large number of particles'
          stop
       endif
    endif

    !charge
    if(nrank == nroot) n0 = np2(nxs,1)
    call MPI_BCAST(n0,1,mnpi,nroot,ncomw,nerr)
    q(1) = fpi*dsqrt(r(1)/(4.0*pi*n0))
    q(2) = -q(1)

    !Magnetic field strength
    b0 = fgi*r(1)*c/q(1)

    call random_gen__init(nrank)
    call init__loading
    call init__set_field
    call fio__param(np,nsp,np2,nxgs,nxge,nxs,nxe,nxs1,nxe1,bcp, &
                    c,q,r,vti,vte,va,rtemp,fpe,fge,             &
                    ldb,delt,delx,dir,file9,                    &
                    nroot,nrank)

  end subroutine init__set_param


  subroutine init__loading

    use boundary, only : boundary__particle

    integer :: i, ii, isp
    real(8) :: sd, sd2, r1, r2

    !particle position
    isp=1
    do i=nxs,nxe+bcp
       do ii=1,np2(i,isp)
          call random_number(r1)
          up(1,ii,i,1) = dble(i)+delx*r1
          up(1,ii,i,2) = up(1,ii,i,1)
       enddo
    enddo

    !velocity
    !Maxwellian distribution
    do isp=1,nsp
       if(isp .eq. 1) then 
          sd = vti/dsqrt(2.0D0)
          sd = sd/sqrt(1.-(sd/c)**2)
          
          do i=nxs,nxe+bcp
             do ii=1,np2(i,isp)
                call random_gen__bm(r1,r2)
                up(2,ii,i,isp) = sd*r1

                call random_gen__bm(r1,r2)
                up(3,ii,i,isp) = sd*r1
                up(4,ii,i,isp) = sd*r2
             enddo
          enddo
       endif
       if(isp .eq. 2) then
          sd = vte/dsqrt(2.0D0)

          do i=nxs,nxe+bcp
             do ii=1,np2(i,isp)
                sd2 = sd/sqrt(1.-(sd/c)**2)
                call random_gen__bm(r1,r2)
                up(2,ii,i,isp) = sd2*r1

                sd2 = sd*sqrt(t_ani)
                sd2 = sd2/sqrt(1.-(sd2/c)**2)
                call random_gen__bm(r1,r2)
                up(3,ii,i,isp) = sd2*r1
                up(4,ii,i,isp) = sd2*r2
             enddo
          enddo
       endif
    enddo

    call boundary__particle(up,                                            &
                            np,nsp,np2,nxgs,nxge,nxs,nxe,nxs1,nxe1,bc,bcp, &
                            nup,ndown,nstat,mnpi,mnpr,ncomw,nerr)

  end subroutine init__loading


  subroutine init__set_field

    use boundary, only : boundary__field

    integer :: i

    !magnetic field
    do i=nxs,nxe+bcp
       uf(1,i) = b0
    enddo
    do i=nxs,nxe
       uf(2,i) = 0.0
       uf(3,i) = 0.0
    enddo

    !electric field
    do i=nxs,nxe
       uf(4,i) = 0.0
    enddo
    do i=nxs,nxe+bcp
       uf(5,i) = 0.0
       uf(6,i) = 0.0
    enddo

    call boundary__field(uf,                       &
                         nxs,nxe,nxs1,nxe1,bc,bcp, &
                         nup,ndown,nroot,nproc,nrank,mnpr,nstat,ncomw,nerr)

  end subroutine init__set_field


end module init
