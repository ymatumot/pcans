module init

  use const
  use mpi_set
  use random_gen

  implicit none

  private

  public :: init__set_param, init__inject

  integer, allocatable, public :: np2(:,:)
  real(8),              public :: q(nsp), r(nsp)
  real(8), allocatable, public :: uf(:,:,:)
  real(8), allocatable, public :: up(:,:,:,:)
  real(8), allocatable, public :: gp(:,:,:,:)
 
  real(8), save :: v0, gam0, b0, vti, vte, delt


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
    v0   = va*ma
    gam0 = 1./sqrt(1.-(v0/c)**2)
    ldmp = ldmp*c/fpi
   
    np2(nys:nye,1) = n0*(nxe+bc-nxs+1)
    np2(nys:nye,2) = np2(nys:nye,1)

    !charge
    q(1) = fpi*sqrt(r(1)/(4.0*pi*n0/delx**2))
    q(2) = -q(1)

    !MAGNETIC FIELD STRENGTH
    b0 = fgi*r(1)*c/q(1)

    !INITIALIZATION OF SUBROUTINES
    call boundary__init(np,nsp,&
                        nxgs,nxge,nygs,nyge,nxs,nxe,nys,nye,bc,    &
                        nup,ndown,mnpi,mnpr,ncomw,nerr,nstat,delx, &
                        u0x=(/v0*gam0,0.0D0/))
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
    call fio__param(np2,0.5*r(1)*vti**2,rtemp,fpe,fge,ldb)

  end subroutine init__set_param


  subroutine init__loading

    integer :: i, j, ii, isp
    real(8) :: sd, r1, r2, gamp, v1, gam1
    !INITIAL VELOCITY PROFILE
    real(8) :: vfunc,x0,eps=1d-40
    vfunc(x0) = -v0*tanh( (x0-nxge*delx)/(ldmp+eps) )

    !MAGNETIC & ELECTRIC FIELDS
    do j=nys-2,nye+2
    do i=nxs-2,nxe+2
       uf(1,i,j) = 0.0D0
       uf(2,i,j) = 0.0D0
       uf(3,i,j) = b0
       uf(4,i,j) = 0.0D0
       uf(5,i,j) = vfunc(i*delx)*b0/c
       uf(6,i,j) = 0.0D0
    enddo
    enddo

    !PARTICLE POSITION
    isp = 1
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
             up(4,ii,j,isp) = sd*r2

             call random_gen__bm(r1,r2)
             up(5,ii,j,isp) = sd*r1

             gamp = sqrt(1.D0+(up(3,ii,j,isp)**2+up(4,ii,j,isp)**2+up(5,ii,j,isp)**2)/(c*c))

             v1 = vfunc(up(1,ii,j,isp))
             gam1 = 1.0D0/dsqrt(1.0D0-v1**2/c**2)

             ! DENSITY FIX: ZENITANI, PHYS. PLASMAS 22, 042116 (2015)
             call random_number(r1)
             if(up(3,ii,j,isp)*v1 >= 0.)then
                up(3,ii,j,isp) = (+up(3,ii,j,isp)+v1*gamp)*gam0
             else
                if(r1 < (-v1*up(3,ii,j,isp)/gamp))then
                   up(3,ii,j,isp) = (-up(3,ii,j,isp)+v1*gamp)*gam0
                else
                   up(3,ii,j,isp) = (+up(3,ii,j,isp)+v1*gamp)*gam0
                endif
             endif
          enddo
       enddo

    enddo

  end subroutine init__loading


  subroutine init__inject

    integer :: isp, ii, ii2, ii3, j, dn
    real(8) :: sd, r1, r2, dx, gamp

    !INJECT PARTICLES IN X=NXS~NXS+V0*DT
    dx  = v0*delt/delx
    dn  = int(n0*dx+0.5)

    do j=nys,nye
       do ii=1,dn
          ii2 = np2(j,1)+ii
          ii3 = np2(j,2)+ii
          call random_number(r1)
          up(1,ii2,j,1) = nxs*delx+r1*dx*delx
          up(1,ii3,j,2) = up(1,ii2,j,1)

          call random_number(r1)
          up(2,ii2,j,1) = dble(j)*delx+delx*r1
          up(2,ii3,j,2) = up(2,ii2,j,1)
       enddo
    enddo

    !VELOCITY
    !MAXWELLIAN DISTRIBUTION
    do isp=1,nsp

       if(isp == 1) then 
          sd = vti/sqrt(2.0D0)
       endif
       if(isp == 2) then
          sd = vte/sqrt(2.0D0)
       endif

       do j=nys,nye
          do ii=np2(j,isp)+1,np2(j,isp)+dn
             call random_gen__bm(r1,r2)
             up(3,ii,j,isp) = sd*r1
             up(4,ii,j,isp) = sd*r2

             call random_gen__bm(r1,r2)
             up(5,ii,j,isp) = sd*r1

             gamp = sqrt(1.D0+(up(3,ii,j,isp)**2+up(4,ii,j,isp)**2+up(5,ii,j,isp)**2)/(c*c))

             call random_number(r1)
             if(up(3,ii,j,isp)*v0 >= 0.)then
                up(3,ii,j,isp) = (+up(3,ii,j,isp)+v0*gamp)*gam0
             else
                if(r1 < (-v0*up(3,ii,j,isp)/gamp))then
                   up(3,ii,j,isp) = (-up(3,ii,j,isp)+v0*gamp)*gam0
                else
                   up(3,ii,j,isp) = (+up(3,ii,j,isp)+v0*gamp)*gam0
                endif
             endif
          enddo
       enddo

    enddo

    do isp=1,nsp
       do j=nys,nye
          np2(j,isp) = np2(j,isp)+dn
       enddo
    enddo

    !SET EX AND BZ
    do j=nys-2,nye+2
       uf(3,nxs-2:nxs+1,j) = b0
       uf(5,nxs-2:nxs+1,j) = v0*b0/c
    enddo

  end subroutine init__inject


end module init
