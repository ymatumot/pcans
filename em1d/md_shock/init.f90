module init

  use const
  use random_gen

  implicit none

  private

  public :: init__set_param, init__inject

  integer, public :: np2(1:nx+bc,nsp)
  real(8), public :: delx, delt
  real(8), public :: c
  real(8), public :: uf(6,0:nx+1)
  real(8), public :: up(4,np,1:nx+bc,nsp)
  real(8), public :: gp(4,np,1:nx+bc,nsp)
  real(8), public :: q(nsp), r(nsp)
  real(8), public :: den(0:nx+1,nsp),vel(0:nx+1,3,nsp),temp(0:nx+1,3,nsp)
  real(8)         :: v0, u0, bx0, by0, bz0, vt0(nsp)


contains


  subroutine init__set_param

    use fio, only : fio__init, fio__input, fio__param
    use boundary, only : boundary__init
    use particle, only : particle__init
    use field, only : field__init
    use mom_calc, only : mom_calc__init

    real(8) :: fpi, fge, fgi, vti, vae, vai, b0
    character(len=64) :: file11

    r(1)  = rmass
    r(2)  = 1.0
    vti   = vte * sqrt(betai/(betae*rmass))
    c     = vte / sqrt(sigma * 0.5*betae)
    fge   = sqrt(sigma)
    fpi   = fpe / sqrt(rmass)
    fgi   = fge / rmass
    vae   = sqrt(sigma) * c
    vai   = sqrt(sigma/rmass) * c
    delx  = vte/fpe/rdbl
    delt  = cfl*delx/c

    ! FOR UPSTREAM BOUNDARY CONDITION
    vt0(1) = vti
    vt0(2) = vte
    v0   = -ma * vai
    u0   = v0/sqrt(1.0D0-(v0/c)**2)

    ! NUMBER OF PARTICLES
    np2(1:nx+bc,1:2) = n0
    if(max(np2(1,1), np2(nx+bc,1), np) > np)then
       write(*,*)'Too large number of particles'
       stop
    endif

    ! MAGNETIC FIELD STRENGTH
    b0   = sqrt(4*pi*r(2)*n0/delx) * vae
    bx0  = b0 * cos(theta)
    by0  = 0.0d0
    bz0  = b0 * sin(theta)

    ! ELEMENTARY CHARGE
    q(1) = fpi * sqrt(r(1)/(4.0*pi*n0/delx))
    q(2) = -q(1)

    !INTIALIZATION OF SUBROUTINES
    call random_gen__init
    call boundary__init(np,nx,nsp,bc)
    call particle__init(np,nx,nsp,bc,q,r,c,delx,delt)
    call field__init(np,nx,nsp,bc,q,c,delx,delt,gfac)
    call fio__init(np,nx,nsp,bc,q,r,c,delx,delt,pi,dir,dir_mom,dir_psd)
    call mom_calc__init(np,nx,nsp,bc,q,r,c,delx,0.5*delt)
    call random_gen__init

    !READING RESTART DATA IF NECESSARY
    if(it0 /= 0)then
       !start from the past calculation
       write(file11,'(i6.6,a)')it0,'_file10.dat'
       call fio__input(up,uf,np2,it0,file11)
       return
    endif

    call init__loading
    call init__set_field

    !OUTPUT PARAMETERS
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

    integer :: i, ii, isp
    real(8) :: r1, r2, uu

    !PARTICLE POSITION
    isp = 1
    do i=1,nx+bc
       do ii=1,np2(i,isp)
          call random_number(r1)
          up(1,ii,i,1) = dble(i)+delx*r1
          up(1,ii,i,2) = up(1,ii,i,1)
       enddo
    enddo

    !VELOCITY
    !MAXWELLIAN DISTRIBUTION
    do isp=1,nsp
       do i=1,nx+bc
          if( dble(i-1) < lbuf ) then
             uu = u0 * dble(i-1)/lbuf
          else
             uu = u0
          end if
          do ii=1,np2(i,isp)
             call random_gen__bm(r1,r2)
             up(2,ii,i,isp) = vt0(isp)*r1 + uu
             up(3,ii,i,isp) = vt0(isp)*r2

             call random_gen__bm(r1,r2)
             up(4,ii,i,isp) = vt0(isp)*r1
          enddo
       enddo
    enddo

  end subroutine init__loading


  subroutine init__set_field

    integer :: i
    real(8) :: uu, vv

    do i=0,nx+1+bc
       uf(1,i) = bx0
    enddo
    do i=0,nx+1
       uf(2,i) = by0
       uf(3,i) = bz0
    enddo

    !ELECTRIC FIELD
    do i=0,nx+1
       uf(4,i) = 0.0
    enddo
    do i=0,nx+1+bc
       if( dble(i-1) < lbuf ) then
          uu = u0 * dble(i-1)/lbuf
          vv = uu/sqrt(1.0D0+(uu/c)**2)
       else
          vv = v0
       end if
       uf(5,i) = +vv*bz0/c
       uf(6,i) = -vv*by0/c
    enddo

  end subroutine init__set_field


  subroutine init__inject

    integer :: isp, i, ii, ii1, ii2, ii3, dn
    real(8) :: rr, r1, r2, ninj, iinj, finj, x0

    i  = nx-1
    x0 = dble(nx)+v0*delt

    ! number of particles to be injected
    call random_number(rr)
    ninj = n0*abs(v0)*delt
    iinj = floor(ninj)
    finj = ninj - iinj
    dn   = int(iinj) + ceiling(finj-rr)

    ! inject particles
    do ii = 1, dn
       ii1 = np2(i,1) + ii
       ii2 = np2(i,2) + ii
       ! position
       call random_number(rr)
       rr = rr*v0*delt + dble(i+1)
       up(1,ii1,i,1) = rr
       up(1,ii2,i,2) = rr

       ! velocity
       do isp = 1, nsp
          ii3 = np2(i,isp) + ii
          call random_gen__bm(r1,r2)
          up(2,ii3,i,isp) = vt0(isp)*r1 + u0
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
    uf(5,i) =+v0*bz0/c
    uf(6,i) =-v0*by0/c

  end subroutine init__inject

end module init
