module init

  use const
  use mpi_set
  use random_gen

  implicit none

  private

  public :: init__set_param

  integer, allocatable, public :: np2(:,:)
  real(8),              public :: q(nsp), r(nsp)
  real(8), allocatable, public :: uf(:,:,:)
  real(8), allocatable, public :: up(:,:,:,:)
  real(8), allocatable, public :: gp(:,:,:,:)
 
  real(8) :: b0, vti, vte, delt


contains

  
  subroutine init__set_param

    use fio, only      : fio__init, fio__input, fio__param
    use boundary, only : boundary__init
    use field, only    : field__init
    use particle, only : particle__init

    real(8) :: fgi, fpi, va, fpe, fge, rgi, rge, ldb
    character(len=64) :: fname_restart

!************** MPI settings  *******************!
    call mpi_set__init(nxgs,nxge,nygs,nyge,nproc)
!*********** End of MPI settings  ***************!
!*********** MEMORY ALLOCATIONS *****************!
    allocate(np2(nys:nye,nsp))
    allocate(uf(6,nxs-2:nxe+2,nys-2:nye+2))
    allocate(up(5,np,nys:nye,nsp))
    allocate(gp(5,np,nys:nye,nsp))
!*********** End of ALLOCATIONS  ***************!

!**** SETTING OTHER NUMERICAL & PHYSICAL CONSTANTS ****!
    r(2) = 1.0D0      ! ELECTRON MASS
    r(1) = r(2)*mr    ! ION MASS
    delt = cfl*delx/c ! TIME STEP SIZE
    ldb  = delx*rdbl
    fpe  = sqrt(beta*rtemp)*c/(sqrt(2.D0)*alpha*ldb)
    fge  = fpe/alpha
    fgi  = fge*r(2)/r(1)
    fpi  = fpe*sqrt(r(2)/r(1))
    va   = fge/fpe*c*sqrt(r(2)/r(1))
    rge  = alpha*ldb*sqrt(2.D0)
    rgi  = rge*sqrt(r(1)/r(2))/sqrt(rtemp)
    vte  = rge*fge
    vti  = vte*sqrt(r(2)/r(1))/sqrt(rtemp)

    np2(nys:nye,1) = n0*(nxe+bc-nxs+1)
    np2(nys:nye,2) = np2(nys:nye,1)

    !ELEMENTARY CHARGE
    q(1) = fpi*sqrt(r(1)/(4.0*pi*n0/delx**2))
    q(2) = -q(1)

    !MAGNETIC FIELD STRENGTH
    b0 = fgi*r(1)*c/q(1)

    !INITIALIZATION OF SUBROUTINES
    call boundary__init(np,nsp,&
                        nxgs,nxge,nygs,nyge,nxs,nxe,nys,nye,bc, &
                        nup,ndown,mnpi,mnpr,ncomw,nerr,nstat,delx)
    call particle__init(np,nsp,&
                        nxs,nxe,nys,nye,nsfo, &
                        q,r,c,delx,delt)
    call field__init(np,nsp,&
                     nxs,nxe,nys,nye,nsfo,bc, &
                     q,c,delx,delt,gfac,      &
                     nup,ndown,mnpr,opsum,ncomw,nerr,nstat)
    call fio__init(np,nsp,&
                   nxgs,nxge,nygs,nyge,nxs,nxe,nys,nye,nsfo,bc, &
                   q,r,c,delx,delt,gfac,n0,                     &
                   fname_param,fname_energy,dir,                &
                   nproc,nrank,nroot,nup,ndown,mnpr,opsum,ncomw,nerr,nstat)
    call random_gen__init(nrank)

    !START FROM RESTART DATA IF AVAILABLE
    if(it0 /= 0)then
       !START FROM THE PAST CALCULATION
       write(fname_restart,'(i6.6,a,i3.3,a)') it0,'_rank=',nrank,'.dat'
       call fio__input(up,uf,np2,it0,fname_restart)
       return
    endif

    !SETUP FOR PARTICLES & FIELDS
    call init__loading
    call init__set_field
    call fio__param(np2,0.5*r(1)*vti**2,rtemp,fpe,fge,ldb)

  end subroutine init__set_param


  subroutine init__loading

    integer :: j, ii, isp
    real(8) :: sd, r1, r2

    !PARTICLE POSITION
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

    !VELOCITY
    !MAXWELLIAN DISTRIBUTION
    do isp=1,nsp
       if(isp .eq. 1) then 
          sd = vti/sqrt(2.0D0)
       endif
       if(isp .eq. 2) then
          sd = vte/sqrt(2.0D0)
       endif

       do j=nys,nye
          do ii=1,np2(j,isp)
             call random_gen__bm(r1,r2)
             up(3,ii,j,isp) = sd*r1

             call random_gen__bm(r1,r2)
             up(4,ii,j,isp) = sd*r1
             up(5,ii,j,isp) = sd*r2
          enddo
       enddo
    enddo

  end subroutine init__loading


  subroutine init__set_field

    integer :: i, j

    !MAGNETIC & ELECTRIC FIELDS
    do j=nys-2,nye+2
    do i=nxs-2,nxe+2
       uf(1,i,j) = 0.0
       uf(2,i,j) = 0.0
       uf(3,i,j) = b0
       uf(4,i,j) = 0.0
       uf(5,i,j) = 0.0
       uf(6,i,j) = 0.0
    enddo
    enddo

  end subroutine init__set_field


end module init
