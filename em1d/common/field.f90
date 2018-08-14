module field

  implicit none

  private

  public :: field__init
  public :: field__fdtd_i

  logical, save              :: is_init = .false.
  integer, save              :: np, nx, nsp, bc
  real(8), save              :: c, delx, delt, gfac, pi, fac1, fac2, ifac1
  real(8), save, allocatable :: q(:), gf(:,:)


contains


  subroutine field__init(npin,nxin,nspin,bcin,qin,cin,delxin,deltin,gfacin)
    integer, intent(in) :: npin, nxin, nspin, bcin
    real(8), intent(in) :: qin(nspin), cin, delxin, deltin, gfacin  

    np   = npin
    nx   = nxin
    nsp  = nspin
    bc   = bcin
    allocate(q(nsp))
    q    = qin
    c    = cin
    delx = delxin
    delt = deltin
    gfac = gfacin
    pi   = 4.0*atan(1.0)
    allocate(gf(6,0:nx+1))
    gf(1:6,0:nx+1) = 0.0
    fac1  = 2.0+(delx/(c*delt*gfac))**2
    ifac1 = 1./fac1
    fac2  = (delx/(c*delt*gfac))**2
    is_init = .true.

  end subroutine field__init

  
  subroutine field__fdtd_i(uf,up,gp,np2)

    use boundary, only : boundary__curre, boundary__field

    !Implicit EM field solver
    integer, intent(in)    :: np2(1:nx+bc,nsp)
    real(8), intent(in)    :: gp(4,np,1:nx+bc,nsp)
    real(8), intent(inout) :: up(4,np,1:nx+bc,nsp)
    real(8), intent(inout) :: uf(6,0:nx+1)
    integer                :: ii, i, isp
    real(8)                :: uj(3,-1:nx+2), gkl(6,0:nx+1)
    real(8)                :: f1, f2, f3, rotb2, rotb3, rote2, rote3

    if(.not.is_init)then
       write(6,*)'Initialize first by calling field__init()'
       stop
    endif

    call ele_cur(uj,up,gp,np2)
    call boundary__curre(uj)

    !calculation
    !< gk =  (c*delt)*rot(b) - (4*pi*delt)*j >
    !< gl =  (c*delt)*rot(e) >
    f1 = c*delt/delx
    f2 = 4.0*pi*delt

    !By,z
    do i=1,nx
       rote3 = -uf(5,i-1)+uf(5,i)
       rote2 = +uf(6,i-1)-uf(6,i)
       gkl(2,i) = -f1*rote2
       gkl(3,i) = -f1*rote3
    enddo
    !Ex
    do i=1,nx
       gkl(4,i) = -f2*uj(1,i)
    enddo
    !Ey,z
    do i=1,nx+bc
       rotb3 = -uf(2,i)+uf(2,i+1)
       rotb2 = +uf(3,i)-uf(3,i+1)
       gkl(5,i) = f1*rotb2-f2*uj(2,i)
       gkl(6,i) = f1*rotb3-f2*uj(3,i)
    enddo

    call boundary__field(gkl)

    f3 = c*delt*gfac/delx
    do i=1,nx
       gkl(2,i) = gkl(2,i)+f3*(-gkl(6,i-1)+gkl(6,i))
       gkl(3,i) = gkl(3,i)-f3*(-gkl(5,i-1)+gkl(5,i))
    enddo

    !solve  < by & bz;  bx = const >
    call cgm(gf,gkl)
    call boundary__field(gf)

    !solve  < ex, ey & ez >
    do i=1,nx
       gf(4,i) = gkl(4,i)   ! if (Ex=0) gf(4,i)=0.
    enddo
    do i=1,nx+bc
       gf(5,i) = gkl(5,i)-f3*(-gf(3,i)+gf(3,i+1))
       gf(6,i) = gkl(6,i)+f3*(-gf(2,i)+gf(2,i+1))
    enddo

    call boundary__field(gf)

    !update fields
    uf(1:6,0:nx+1) = uf(1:6,0:nx+1)+gf(1:6,0:nx+1)

    !update particle
    do isp=1,nsp
       do i=1,nx+bc
          do ii=1,np2(i,isp)
             up(1,ii,i,isp) = gp(1,ii,i,isp)
             up(2,ii,i,isp) = gp(2,ii,i,isp)
             up(3,ii,i,isp) = gp(3,ii,i,isp)
             up(4,ii,i,isp) = gp(4,ii,i,isp)
          enddo
       enddo
    enddo


  end subroutine field__fdtd_i


  subroutine ele_cur(uj,up,gp,np2)

    integer, intent(in)  :: np2(1:nx+bc,nsp)
    real(8), intent(in)  :: up(4,np,1:nx+bc,nsp), gp(4,np,1:nx+bc,nsp)
    real(8), intent(out) :: uj(3,-1:nx+2)
    integer :: ii, i, isp, ih, i1 ,i2
    real(8) :: x2, xh, xr, qvx1, qvx2, idelt, idelx, gam
    real(8) :: dx, dxm

    !memory clear
    uj(1:3,-1:nx+2) = 0.0D0

    !------ Charge Conservation Method for Jx ---------!
    idelt = 1.D0/delt
    idelx = 1.D0/delx
    do isp=1,nsp
       do i=1,nx+bc
          do ii=1,np2(i,isp)

             x2 = gp(1,ii,i,isp)

             !reflective boundary condition in x
             if(bc == -1)then
                if(x2 < 1)then
                   x2  = 2.-x2
                endif
                if(x2 > nx)then
                   x2  = 2.*nx-x2
                endif
             endif

             i1  = floor(up(1,ii,i,isp)*idelx-0.5)
             i2  = floor(x2*idelx-0.5)

             xh  = 0.5*(up(1,ii,i,isp)+x2)

             if(i1==i2)then
                xr = xh
             else
                xr = (max(i1,i2)+0.5)*delx
             endif

             qvx1 = q(isp)*(xr-up(1,ii,i,isp))*idelt
             qvx2 = q(isp)*(x2-xr)*idelt

             !Jx
             uj(1,i1+1) = uj(1,i1+1)+qvx1
             uj(1,i2+1) = uj(1,i2+1)+qvx2

             !Jy and Jz
             gam = 1./sqrt(1.0D0+(+gp(2,ii,i,isp)*gp(2,ii,i,isp) &
                                  +gp(3,ii,i,isp)*gp(3,ii,i,isp) &
                                  +gp(4,ii,i,isp)*gp(4,ii,i,isp) &
                                 )/(c*c))

             ih = floor(xh-0.5)
             dx = xh-0.5-ih
             dxm = 1.-dx

             uj(2,ih  ) = uj(2,ih  )+q(isp)*gp(3,ii,i,isp)*gam*dxm
             uj(3,ih  ) = uj(3,ih  )+q(isp)*gp(4,ii,i,isp)*gam*dxm
             uj(2,ih+1) = uj(2,ih+1)+q(isp)*gp(3,ii,i,isp)*gam*dx 
             uj(3,ih+1) = uj(3,ih+1)+q(isp)*gp(4,ii,i,isp)*gam*dx 
          enddo
       enddo
    enddo

    if(bc == 0)then
       do i=-1,nx+2
          uj(1:3,i) = uj(1:3,i)*idelx
       enddo
    else if(bc == -1)then
       i=1
       uj(1,i) = uj(1,i)*2.*idelx
       do i=2,nx-1
          uj(1,i) = uj(1,i)*idelx
       enddo
       i=nx
       uj(1,i) = uj(1,i)*2.*idelx

       do i=0,nx
          uj(2,i) = uj(2,i)*idelx
          uj(3,i) = uj(3,i)*idelx
       enddo
    else
       write(*,*)'choose bc=0 (periodic) or bc=-1 (reflective)'
       stop
    endif

  end subroutine ele_cur


  subroutine cgm(gb,gl)

    !-----------------------------------------------------------------------
    !  #  conjugate gradient method 
    !-----------------------------------------------------------------------

    use boundary, only : boundary__phi

    real(8), intent(in)    :: gl(6,0:nx+1)
    real(8), intent(inout) :: gb(6,0:nx+1)
    integer, parameter     :: ite_max = 1000 ! maximum number of iteration
    integer :: i, l, ite
    real(8), parameter     :: err = 1d-6
    real(8) :: eps, sumr, sum, sum1, sum2, av, bv
    real(8) :: x(0:nx+1), b(0:nx+1), r(0:nx+1), p(0:nx+1), ap(0:nx+1)


    do l=2,3

       ! initial guess
       ite = 0
       sum = 0.0
       do i=1,nx
          x(i) = gb(l,i)
          b(i) = fac2*gl(l,i)
          sum = sum+b(i)*b(i)
       enddo

       eps = sqrt(sum)*err
       
       call boundary__phi(x)

       sumr = 0.0
       do i=1,nx
          r(i) = b(i)+x(i-1)-fac1*x(i)+x(i+1)
          p(i) = r(i)
          sumr = sumr+r(i)*r(i)
       enddo

       if(sqrt(sumr) > eps)then
       
          do while(sum > eps)
             
             ite = ite+1

             call boundary__phi(p)
       
             do i=1,nx
                ap(i) = -p(i-1)+fac1*p(i)-p(i+1)
             enddo
             sumr = 0.0
             sum2 = 0.0
             do i=1,nx
                sumr = sumr+r(i)*r(i)
                sum2 = sum2+p(i)*ap(i)
             enddo
             av = sumr/sum2
             
             x(1:nx) = x(1:nx)+av*p(1:nx)
             r(1:nx) = r(1:nx)-av*ap(1:nx)
             
             sum = sqrt(sumr)
             if(ite >= ite_max) then
                write(6,*)'********** stop at cgm after ite_max **********'
                stop
             endif
             
             sum1 = 0.0
             do i=1,nx
                sum1 = sum1+r(i)*r(i)
             end do
             bv = sum1/sumr
             
             p(1:nx) = r(1:nx)+bv*p(1:nx)

          enddo

       endif

       gb(l,1:nx) = x(1:nx)

    end do
    
  end subroutine cgm


end module field
