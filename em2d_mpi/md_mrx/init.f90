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
 
  real(8) :: x0, y0, vti, vte, b0, delt


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

! HERE WE USE THE PRESSURE BALANCE,
!   B**2 / 8*pi = n [ kT_e + kT_i ] = n kT_e * ((1+rtemp)/rtemp)
!               = 0.5*me*n*(vte**2) * ((1+rtemp)/rtemp)
!   (vte**2) * ((1+rtemp)/rtemp) = B**2 / 4*pi*me*n = (c/alpha)**2
    vte = sqrt(rtemp)*c/(sqrt(1+rtemp)*alpha)
    vti = vte*sqrt(r(2)/r(1))/sqrt(rtemp)
! WE USE THE FOLLOWING RELATION TO OBTAIN FPE,
!   ldb = sqrt( kT_e / 4 pi n e**2 )
!       = sqrt( 0.5*me*(vte**2) / 4 pi n e**2 )
!       = vte/fpe/sqrt(2)
    fpe = vte/ldb/sqrt(2.d0)
    fpi = fpe*sqrt(r(2)/r(1))
    fge = fpe/alpha
    fgi = fge*r(2)/r(1)
    va  = fge/fpe*c*sqrt(r(2)/r(1))
    rge = vte/fge
    rgi = vti/fgi

    ! ELEMENTARY CHARGE
    q(1) = fpi*sqrt(r(1)/(4.d0*pi*ncs))
    q(2) = -q(1)  ! -fpe*sqrt(r(2)/(4.d0*pi*ncs))

    ! MAGNETIC FIELD STRENGTH
    b0  = fgi*r(1)*c/q(1)
    
    ! POSITION OF THE X-POINT
    x0  = 0.5*(nxge+nxgs)*delx 
    y0  = 0.5*(nyge-nygs)*delx
    ! CURRENT SHEET THICKNESS
    lcs = lcs * c/fpi

    ! NUMBER OF PARTICLES IN EACH CELL IN Y
    np2(nys:nye,1:nsp) = ceiling( nbg*(nxge-nxgs)*delx**2 + ncs*2*lcs*delx + 1d-6)

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
                   q,r,c,delx,delt,gfac,ncs,                    &
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
    call fio__param(np2,0.5*r(1)*vti**2,rtemp,fpe,fge,ldb)

  end subroutine init__set_param


  subroutine init__loading

    integer            :: i, j, ii, ibg
    real(8), parameter :: e1 = 0.12d0  ! B1/B0: I recommend ~< O(0.5*nbg/ncs).
    real(8)            :: r1, r2
    real(8)            :: jz, density
    real(8)            :: b_harris, bx_pert, by_pert
    real(8)            :: x, y
    real(8)            :: f1, f2
    real(8)            :: sdi, sde  ! for Box-Mullter method

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
    do j=nys-2,nye+2
    do i=nxs-2,nxe+2
       uf(1,i,j) = 0.0D0
       uf(1,i,j) = uf(1,i,j) + bx_pert(i*delx,j*delx)
       uf(2,i,j) = b_harris(i*delx)
       uf(2,i,j) = uf(2,i,j) + by_pert(i*delx,j*delx)
       uf(3,i,j) = 0.d0
       uf(4,i,j) = 0.d0
       uf(5,i,j) = 0.d0
       uf(6,i,j) = 0.d0
    enddo
    enddo
    ! ---------------- Electromagnetic field ------------------

    ! ---------------- Particles ------------------
    f1 =  1.d0 / ( (1.d0+rtemp) * q(1) )
    f2 = rtemp / ( (1.d0+rtemp) * q(2) )
    ! a factor of 1/sqrt(2) for Box-Muller method
    sdi = vti/sqrt(2.)
    sde = vte/sqrt(2.)

    do j=nys,nye

       ibg = floor( nbg*(nxe+bc-nxs+1)*delx**2 + 1d-6)
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
          up(3,ii,j,1) = sdi*r1
          up(4,ii,j,1) = sdi*r2
          call random_gen__bm(r1,r2)
          up(5,ii,j,1) = sdi*r1 + f1*jz(up(1,ii,j,1),up(2,ii,j,1))/density(up(1,ii,j,1))
          up(3,ii,j,2) = sde*r2
          call random_gen__bm(r1,r2)
          up(4,ii,j,2) = sde*r1
          up(5,ii,j,2) = sde*r2 + f2*jz(up(1,ii,j,2),up(2,ii,j,2))/density(up(1,ii,j,2))
       enddo

    enddo

  end subroutine init__loading


end module init
