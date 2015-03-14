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
  real(8)                   :: pi

  integer :: n0
  real(8) :: b0, vt0(nsp), nbeam, vbeam, tbeam

contains


  subroutine init__set_param

    use fio, only : fio__input, fio__param

    real(8) :: rmass, fpe, fpi, fge, fgi, vte, vti, vae, vai, betae, betai
    real(8) :: sigma, Ma, theta
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
    itmax  = 400000
    intvl1 = 1000
    intvl2 = 1000
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

    pi   = 4.0*atan(1.0)

    !
    ! n0    : number of particles / cell
    ! rmass : mass ratio
    ! sigma : (wce/wpe)^2
    ! betae : electron beta
    ! betai : ion beta
    ! nbeam : density ratio (beam / total)
    ! vbeam : relative streaming velocity
    ! tbeam : temperature ratio (beam / core)
    !
    n0    = 250
    rmass = 25.0
    sigma = 0.04
    betae = 0.5
    betai = 0.03125

    vte   = 1.0
    vti   = vte * dsqrt(betai/(betae*rmass))
    c     = vte / dsqrt(sigma * 0.5*betae)
    fpe   = 1.0
    fge   = dsqrt(sigma)
    fpi   = fpe / dsqrt(rmass)
    fgi   = fge / rmass
    vae   = dsqrt(sigma) * c
    vai   = dsqrt(sigma/rmass) * c
    b0    = dsqrt(4*pi*n0) * vae
    delx  = vte / fpe
    delt  = 0.5*delx/c

    ! beam parameters
    nbeam = 0.2
    vbeam = 3.0 * vai
    tbeam = 1.0
    vt0(1) = vti
    vt0(2) = vte

    np2(1:nx+bc,1:2) = n0

    if(max(np2(1,1), np2(nx+bc,1), np) > np)then
       write(*,*)'Too large number of particles'
       stop
    endif

    ! mass and charge
    r(1) = rmass
    r(2) = 1.0
    q(1) = fpi * dsqrt(r(1)/(4.0*pi*n0))
    q(2) = -q(1)

    call random_gen__init
    call init__loading
    call init__set_field

    ! output parameters
    open(9,file=trim(dir)//trim(file9),status='unknown')

    write(9,"(A30, 2x, i10)") "number of grids : ", nx
    write(9,"(A30, 2x, i10)") "number of particles/cell : ", n0
    write(9,"(A30, 2x, es10.3)") "speed of light : ", c
    write(9,"(A30, 2x, es10.3)") "time step : ", delt
    write(9,"(A30, 2x, es10.3)") "grid size : ", delx
    write(9,"(A30, 2x, es10.3)") "beam density ratio : ", nbeam
    write(9,"(A30, 2x, es10.3)") "relative velocity : ", vbeam
    write(9,"(A30, 2x, es10.3)") "beam temperature ratio : ", tbeam
    write(9,"(A30, 2x, es10.3)") "Debye length : ", vte/fpe
    write(9,"(A30, 2x, es10.3)") "ion mass : ", r(1)
    write(9,"(A30, 2x, es10.3)") "electron mass : ", r(2)
    write(9,"(A30, 2x, es10.3)") "ion charge : ", q(1)
    write(9,"(A30, 2x, es10.3)") "electron charge : ", q(2)
    write(9,"(A30, 2x, es10.3)") "ion plasma freq : ", fpi
    write(9,"(A30, 2x, es10.3)") "electron plasma freq : ", fpe
    write(9,"(A30, 2x, es10.3)") "ion gyro freq : ", fgi
    write(9,"(A30, 2x, es10.3)") "electron gyro freq : ", fge
    write(9,"(A30, 2x, es10.3)") "ion thermal velocity : ", vti
    write(9,"(A30, 2x, es10.3)") "electron thermal velocity : ", vte
    write(9,"(A30, 2x, es10.3)") "ion Alfven velocity : ", vai
    write(9,"(A30, 2x, es10.3)") "electron Alfven velocity : ", vae

    close(9)

  end subroutine init__set_param


  subroutine init__loading

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
       nb  = int( np2(i,isp) * nbeam )
       vb1 =+vbeam * (1 - nbeam)
       vb2 =-vbeam * nbeam
       vt1 = dsqrt(2.0 * tbeam / (1.0 + tbeam)) * vt0(isp)
       vt2 = dsqrt(2.0         / (1.0 + tbeam)) * vt0(isp)
       ! beam componennt
       do ii=1,nb
          call random_gen__bm(r1,r2)
          up(2,ii,i,isp) = vt1*r1 + vb1

          call random_gen__bm(r1,r2)
          up(3,ii,i,isp) = vt1*r1
          up(4,ii,i,isp) = vt1*r2
       end do
       ! core component
       do ii=nb+1,np2(i,isp)
          call random_gen__bm(r1,r2)
          up(2,ii,i,isp) = vt2*r1 + vb2

          call random_gen__bm(r1,r2)
          up(3,ii,i,isp) = vt2*r1
          up(4,ii,i,isp) = vt2*r2
       end do
    end do

    ! electrons
    isp = 2
    do i=1,nx+bc
       do ii=1,np2(i,isp)
          call random_gen__bm(r1,r2)
          up(2,ii,i,isp) = vt0(isp)*r1

          call random_gen__bm(r1,r2)
          up(3,ii,i,isp) = vt0(isp)*r1
          up(4,ii,i,isp) = vt0(isp)*r2
       end do
    end do

  end subroutine init__loading


  subroutine init__set_field

    integer :: i

    !magnetic field
    do i=0,nx+1+bc
       uf(1,i) = b0
    enddo
    do i=0,nx+1
       uf(2,i) = 0.0
       uf(3,i) = 0.0
    enddo

    !electric field
    do i=0,nx+1
       uf(4,i) = 0.0
    enddo
    do i=0,nx+1+bc
       uf(5,i) = 0.0
       uf(6,i) = 0.0
    enddo

  end subroutine init__set_field


end module init
