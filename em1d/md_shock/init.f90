module init

  use const
  use random_gen

  implicit none

  private

  public :: init__set_param, init__inject

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
  real(8) :: v0, u0, b0, bx0, by0, bz0, vt0(nsp)

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
!   gfac    : implicit factor
!             gfac < 0.5 : unstable
!             gfac = 0.5 : no implicit
!             gfac = 1.0 : full implicit
!*********************************************************************
    itmax  = 200000
    intvl1 = 500
    intvl2 = 500
    dir    = './dat/'
    file9  = 'init_param.dat'
    file10 = 'file10.dat'
    file12 = 'energy.dat'
    gfac   = 0.505
    it0    = 0

    pi   = 4.0*atan(1.0)

    !
    ! n0    : number of particles / cell
    ! rmass : mass ratio
    ! sigma : (wce/wpe)^2
    ! betae : electron beta
    ! betai : ion beta
    ! Ma    : Mach number of the flow (in the simulation frame)
    ! theta : upstream magnetic field angle normal to the x axis
    !
    n0    = 64
    rmass = 25.0
    sigma = 0.04
    betae = 1.0/8.0
    betai = 1.0/8.0
    Ma    = 3.0
    theta = 90.0 * pi/180.0

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

    ! for upstream boundary condition
    vt0(1) = vti
    vt0(2) = vte
    v0   = Ma * vai
    u0   = v0/dsqrt(1.0-(v0/c)**2)
    bx0  = b0 * dcos(theta)
    by0  = 0.0
    bz0  = b0 * dsin(theta)

    ! number of particles
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

    if(it0 /= 0)then
       !start from the past calculation
       file11 = '020000_test10.dat'
       call fio__input(up,uf,np2,c,q,r,delt,delx,it0,np,nx,nsp,bc,dir,file11)
       return
    endif

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
    write(9,"(A30, 2x, es10.3)") "Alfven Mach number : ", Ma
    write(9,"(A30, 2x, es10.3)") "Shock angle : ", theta / pi * 180.0
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

    use boundary, only : boundary__particle

    integer :: i, ii, isp
    real(8) :: sd, r1, r2

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
    do isp=1,nsp
       do i=1,nx+bc
          do ii=1,np2(i,isp)
             call random_gen__bm(r1,r2)
             up(2,ii,i,isp) = vt0(isp)*r1 - u0
             up(3,ii,i,isp) = vt0(isp)*r2

             call random_gen__bm(r1,r2)
             up(4,ii,i,isp) = vt0(isp)*r1
          enddo
       enddo
    enddo

    call boundary__particle(up,np2)

  end subroutine init__loading


  subroutine init__set_field

    use boundary, only : boundary__field

    integer :: i

    !magnetic field
    do i=1,nx+bc
       uf(1,i) = bx0
    enddo
    do i=1,nx
       uf(2,i) = by0
       uf(3,i) = bz0
    enddo

    !electric field
    do i=1,nx
       uf(4,i) = 0.0
    enddo
    do i=1,nx+bc
       uf(5,i) =-v0*bz0/c
       uf(6,i) =+v0*by0/c
    enddo

    call boundary__field(uf)

  end subroutine init__set_field


  subroutine init__inject
    integer :: isp, i, ii, ii1, ii2, ii3, dn
    real(8) :: rr, r1, r2, ninj, iinj, finj, x0

    i  = nx-1
    x0 =-v0*delt + dble(nx)

    ! number of particles to be injected
    call random_number(rr)
    ninj = n0*v0*delt
    iinj = floor(ninj)
    finj = ninj - iinj
    dn   = int(iinj) + ceiling(finj-rr)

    ! inject particles
    do ii = 1, dn
       ii1 = np2(i,1) + ii
       ii2 = np2(i,2) + ii
       ! position
       call random_number(rr)
       rr =-rr*v0*delt + dble(i+1)
       up(1,ii1,i,1) = rr
       up(1,ii2,i,2) = rr

       ! velocity
       do isp = 1, nsp
          ii3 = np2(i,isp) + ii
          call random_gen__bm(r1,r2)
          up(2,ii3,i,isp) = vt0(isp)*r1 - u0
          up(3,ii3,i,isp) = vt0(isp)*r2

          call random_gen__bm(r1,r2)
          up(4,ii3,i,isp) = vt0(isp)*r1

          ! folding back
          if( delt*up(2,ii3,i,isp) > up(1,ii3,i,isp)-dble(i+1) ) then
             up(1,ii3,i,isp) = 2*x0 - up(1,ii3,i,isp)
             up(2,ii3,i,isp) = 2*v0 - up(2,ii3,i,isp)
          end if
       end do
    end do

    ! increase number of particles
    do isp = 1, nsp
       np2(i,isp) = np2(i,isp) + dn
    enddo

    ! ex, by, bz
    i = nx
    uf(2,i) = by0
    uf(3,i) = bz0
    uf(4,i) = 0.0
    ! ey, ez, bx
    i = nx-1
    uf(1,i) = bx0
    uf(5,i) =-v0*bz0/c
    uf(6,i) =+v0*by0/c

  end subroutine init__inject

end module init
