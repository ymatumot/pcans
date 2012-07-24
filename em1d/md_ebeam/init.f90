module init

  use const
  use random_gen

  implicit none

  private

  public :: init__set_param

  integer, public :: np2(1:nx+bc,nsp)
  integer, public :: itmax, it0, intvl1, intvl2, intvl3
  real(8), public :: delx, delt, gfac
  real(8), public :: c
  real(8), public :: uf(6,0:nx+1)
  real(8), public :: up(4,np,1:nx+bc,nsp)
  real(8), public :: q(nsp), r(nsp)
  !gx, gv, are temporal spaces used for the time integration
  real(8), public :: gp(4,np,1:nx,nsp) !just for initialization
  character(len=64), public :: dir
  character(len=64), public :: file10
  character(len=64), public :: file12
  character(len=64), public :: file13, file14
  real(8)                   :: pi, vti, vte, va, rtemp, fpe, fge, rgi, rge, ldb, b0

  real(8) :: nbeam, vbeam, tbeam ! for beam setup

contains


  subroutine init__set_param

    use fio, only : fio__input, fio__param

    real(8) :: fgi, fpi, rmass, betae, betai, sigma, vae, vai
    character(len=64) :: file9 
    character(len=64) :: file11

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
!           : 13~14 - for w-k diagram
!   gfac    : implicit factor
!             gfac < 0.5 : unstable
!             gfac = 0.5 : no implicit
!             gfac = 1.0 : full implicit
!*********************************************************************
    itmax  = 80000
    intvl1 = 500
    intvl2 = 500
    dir    = './dat/'
    file9  = 'init_param.dat'
    file10 = 'file10.dat'
    file12 = 'energy.dat'
    gfac   = 0.505

    it0    = 0
    if(it0 /= 0)then
       !start from the past calculation
       file11 = '002048_file10.dat'
       call fio__input(up,uf,np2,c,q,r,delt,delx,it0,np,nx,nsp,bc,dir,file11)
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
!*********************************************************************
    pi   = 4.0*atan(1.0)
    delx = 1.0
    c    = 200.0
    delt = 0.0025

    !
    ! rmass : mass ratio
    ! sigma : (wce/wpe)^2
    ! betae : electron beta
    ! betai : ion beta
    ! nbeam : density ratio (beam / total)
    ! vbeam : relative streaming velocity
    ! tbeam : temperature ratio (beam / core)
    !
    rmass = 100.0
    sigma = 0.01
    vae   = dsqrt(sigma) * c
    vai   = dsqrt(sigma/rmass) * c
    betae = 0.005
    betai = 0.5
    nbeam = 0.5
    vbeam = 20.0
    tbeam = 1.0

    fpe   = 1.0
    fge   = dsqrt(sigma)
    fpi   = fpe / dsqrt(rmass)
    fgi   = fge / rmass
    va    = vai
    vte   = dsqrt(0.5*betae) * vae
    vti   = dsqrt(0.5*betai) * vai
    rtemp = (vte**2) / (rmass*vti**2)
    ldb   = vte / fpe

    np2(1:nx+bc,1) = 500
    np2(1:nx+bc,2) = np2(1:nx+bc,1)

    if(max(np2(1,1), np2(nx+bc,1), np) > np)then
       write(*,*)'Too large number of particles'
       stop
    endif

    !charge
    r(1) = rmass
    r(2) = 1.0
    q(1) = fpi*dsqrt(r(1)/(4.0*pi*np2(1,1)))
    q(2) = -q(1)

    !Magnetic field strength
    b0 = fgi*r(1)*c/q(1)

    call random_gen__init
    call init__loading
    call init__set_field
    call fio__param(np,nx,nsp,np2,c,q,r,vti,vte,va,rtemp,fpe,fge,ldb,delt,delx,bc,dir,file9)

  end subroutine init__set_param


  subroutine init__loading

    use boundary, only : boundary__particle

    integer :: i, ii, isp, nb
    real(8) :: r1, r2, vb1, vb2, vt1, vt2

    !particle position
    isp = 1
    do i=1,nx+bc
       do ii=1,np2(i,isp)
          call random_number(r1)
          up(1,ii,i,1) = dble(i)+delx*r1
          up(1,ii,i,2) = up(1,ii,i,1)
       enddo
    enddo

    !velocity
    !Maxwellian distribution

    ! ions
    isp = 1
    do i=1,nx+bc
       do ii=1,np2(i,isp)
          call random_gen__bm(r1,r2)
          up(2,ii,i,isp) = vti*r1

          call random_gen__bm(r1,r2)
          up(3,ii,i,isp) = vti*r1
          up(4,ii,i,isp) = vti*r2
       end do
    end do

    ! electrons
    isp = 2
    do i=1,nx+bc
       nb  = int( np2(i,isp) * nbeam )
       vb1 =+vbeam * (1 - nbeam)
       vb2 =-vbeam * nbeam
       vt1 = dsqrt(2.0 * tbeam / (1.0 + tbeam)) * vte
       vt2 = dsqrt(2.0         / (1.0 + tbeam)) * vte
       ! component 1
       do ii=1,nb
          call random_gen__bm(r1,r2)
          up(2,ii,i,isp) = vt1*r1 + vb1

          call random_gen__bm(r1,r2)
          up(3,ii,i,isp) = vt1*r1
          up(4,ii,i,isp) = vt1*r2
       end do
       ! component 2
       do ii=nb+1,np2(i,isp)
          call random_gen__bm(r1,r2)
          up(2,ii,i,isp) = vt2*r1 + vb2

          call random_gen__bm(r1,r2)
          up(3,ii,i,isp) = vt2*r1
          up(4,ii,i,isp) = vt2*r2
       end do
    end do

  end subroutine init__loading


  subroutine init__set_field

    use boundary, only : boundary__field

    integer :: i

    !magnetic field
    do i=1,nx+bc
       uf(1,i) = b0
    enddo
    do i=1,nx
       uf(2,i) = 0.0
       uf(3,i) = 0.0
    enddo

    !electric field
    do i=1,nx
       uf(4,i) = 0.0
    enddo
    do i=1,nx+bc
       uf(5,i) = 0.0
       uf(6,i) = 0.0
    enddo

    call boundary__field(uf)

  end subroutine init__set_field


end module init
