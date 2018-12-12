module init

  use const
  use random_gen

  implicit none

  private

  public :: init__set_param

  integer, public :: np2(1:nx+bc,nsp)
  real(8), public :: delx, delt
  real(8), public :: c
  real(8), public :: uf(6,0:nx+1)
  real(8), public :: up(4,np,1:nx+bc,nsp)
  real(8), public :: q(nsp), r(nsp)
  real(8), public :: gp(4,np,1:nx+bc,nsp)
  real(8), public :: den(0:nx+1,nsp),vel(0:nx+1,3,nsp),temp(0:nx+1,3,nsp)
  real(8)         :: b0, vt0(nsp)


contains


  subroutine init__set_param

    use fio, only : fio__init, fio__input, fio__param
    use boundary, only : boundary__init
    use particle, only : particle__init
    use field, only : field__init
    use mom_calc, only : mom_calc__init

    real(8) :: fpi, fge, fgi, vti, vae, vai
    character(len=64) :: file11

    r(1)  = rmass
    r(2)  = 1.0
    vti   = vte * sqrt(betai/(betae*rmass))
    c     = vte / sqrt(sigma * 0.5*betae)
    fge   = sqrt(sigma)*fpe
    fpi   = fpe / sqrt(rmass)
    fgi   = fge / rmass
    vae   = sqrt(sigma) * c
    vai   = sqrt(sigma/rmass) * c
    delx  = vte/fpe/rdbl
    delt  = cfl*delx/c

    ! BEAM PARAMETERS
    vt0(1) = vti
    vt0(2) = vte

    np2(1:nx+bc,1:2) = n0

    if(max(np2(1,1), np2(nx+bc,1), np) > np)then
       write(*,*)'Too large number of particles'
       stop
    endif

    ! MAGNETIC FIELD STRENGTH
    b0   = sqrt(4.0*pi*r(2)*n0/delx) * vae

    ! ELEMENTARY CHARGE
    q(1) = fpi * sqrt(r(1)/(4.0*pi*n0/delx))
    q(2) = -q(1)

    !INTIALIZATION OF SUBROUTINES
    call random_gen__init
    call boundary__init(np,nx,nsp,bc,delx)
    call particle__init(np,nx,nsp,bc,q,r,c,delx,delt)
    call field__init(np,nx,nsp,bc,q,c,delx,delt,gfac)
    call fio__init(np,nx,nsp,bc,q,r,c,delx,delt,pi,dir,dir_mom,dir_psd)
    call mom_calc__init(np,nx,nsp,bc,q,r,c,delx,0.5*delt)

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
          up(1,ii,i,1) = (dble(i)+r1)*delx
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
       vt1 = sqrt(2.0D0 * tbeam / (1.0D0 + tbeam)) * vt0(isp)
       vt2 = sqrt(2.0D0         / (1.0D0 + tbeam)) * vt0(isp)
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
