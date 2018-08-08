module init

  use const
  use random_gen

  implicit none

  private

  public :: init__set_param

  integer, public :: np2(1:nx+bc,nsp)
  real(8), public :: uf(6,0:nx+1)
  real(8), public :: up(4,np,1:nx+bc,nsp)
  real(8), public :: gp(4,np,1:nx+bc,nsp)
  real(8), public :: q(nsp), r(nsp)
  real(8), public :: delt
  real(8), public :: den(0:nx+1,nsp),vel(0:nx+1,3,nsp),temp(0:nx+1,3,nsp)
  real(8)         :: b0, vte, vti


contains

  
  subroutine init__set_param

    use fio, only : fio__init, fio__input, fio__param
    use boundary, only : boundary__init
    use particle, only : particle__init
    use field, only : field__init
    use mom_calc, only : mom_calc__init

    real(8)           :: va, fpe, fge, rgi, rge, ldb, fgi, fpi
    character(len=64) :: file11

!**** SETTING OTHER NUMERICAL & PHYSICAL CONSTANTS ****!
    r(2) = 1.0D0      ! ELECTRON MASS
    r(1) = r(2)*mr    ! ION MASS
    delt = cfl*delx/c ! TIME STEP SIZE
    ldb  = delx*rdbl
    fpe  = dsqrt(beta*rtemp)*c/(dsqrt(2.D0)*alpha*ldb)
    fge  = fpe/alpha
    fgi  = fge*r(2)/r(1)
    fpi  = fpe*dsqrt(r(2)/r(1))
    va   = fge/fpe*c*dsqrt(r(2)/r(1))
    rge  = alpha*ldb*dsqrt(2.D0)
    rgi  = rge*dsqrt(r(1)/r(2))/dsqrt(rtemp)
    vte  = rge*fge
    vti  = vte*dsqrt(r(2)/r(1))/dsqrt(rtemp)

    np2(1:nx+bc,1) = n0
    np2(1:nx+bc,2) = np2(1:nx+bc,1)
    
    if(max(np2(1,1), np2(nx+bc,1), np) > np)then
       write(*,*)'Too large number of particles'
       stop
    endif

    !CHARGE
    q(1) = fpi*dsqrt(r(1)/(4.0*pi*np2(1,1)))
    q(2) = -q(1)

    !MAGNETIC FIELD STRENGTH
    b0 = fgi*r(1)*c/q(1)

    !INTIALIZATION OF SUBROUTINES
    call random_gen__init
    call boundary__init(np,nx,nsp,bc)
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
    call fio__param(np2,vti,vte,va,rtemp,fpe,fge,ldb,file9)

  end subroutine init__set_param


  subroutine init__loading

    integer :: i, ii, isp
    real(8) :: sd, r1, r2, u0

    u0 = v0/dsqrt(1.-(v0/c)**2)

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
       if(isp == 1) then 
          sd = vti/dsqrt(2.0D0)
       endif
       if(isp == 2) then
          sd = vte/dsqrt(2.0D0)
       endif

       do i=1,nx+bc
          do ii=1,np2(i,isp)
             call random_gen__bm(r1,r2)
             up(2,ii,i,isp) = sd*r1+u0

             call random_gen__bm(r1,r2)
             up(3,ii,i,isp) = sd*r1
             up(4,ii,i,isp) = sd*r2
          enddo
       enddo
    enddo

  end subroutine init__loading


  subroutine init__set_field

    integer :: i

    !UF(1:3): MAGNETIC FIELD
    !UF(4:6): ELECTRIC FIELD
    do i=0,nx+1
       uf(1,i) = b0
       uf(2,i) = 0.0D0
       uf(3,i) = 0.0D0
       uf(4,i) = 0.0D0
       uf(5,i) = 0.0D0
       uf(6,i) = 0.0D0
    enddo

  end subroutine init__set_field


end module init
