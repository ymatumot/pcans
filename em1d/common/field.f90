module field

  implicit none

  private

  public :: field__fdtd_i


contains

  
  subroutine field__fdtd_i(uf,up,gp,np,nx,nsp,np2,bc,q,c,delx,delt,gfac)

    use boundary, only : boundary__curre, boundary__field

    !Implicit EM field solver
    integer, intent(in)    :: np, nx, nsp, bc
    integer, intent(in)    :: np2(1:nx+bc,nsp)
    real(8), intent(in)    :: q(nsp), c, delx, delt, gfac
    real(8), intent(in)    :: gp(4,np,1:nx+bc,nsp)
    real(8), intent(inout) :: up(4,np,1:nx+bc,nsp)
    real(8), intent(inout) :: uf(6,0:nx+1)
    logical, save              :: lflag=.true.
    integer                    :: ii, i, isp
    real(8)                    :: uj(3,-1:nx+2), gkl(6,0:nx+1)
    real(8), save, allocatable :: gf(:,:)
    real(8)                    :: pi, f1, f2, f3, rotb2, rotb3, rote2, rote3

    pi = 4.0*atan(1.0)

    if(lflag)then
       allocate(gf(6,0:nx+1))
       gf(1:3,0:nx+1) = 0.0
       lflag = .false.
    endif
       
!!$    !position at n+1/2
!!$    do isp=1,nsp
!!$       do i=1,nx+bc
!!$          do ii=1,np2(i,isp)
!!$             up(1,ii,i,isp) = 0.5*(up(1,ii,i,isp)+gp(1,ii,i,isp))
!!$             up(2:4,ii,i,isp) = gp(2:4,ii,i,isp)
!!$          enddo
!!$       enddo
!!$    enddo
!!$    call ele_cur2(uj,up,np,nx,nsp,np2,bc,q,c,delx)
    call ele_cur(uj,up,gp,np,nx,nsp,np2,bc,q,c,delx,delt)
    call boundary__curre(uj,nx,bc)

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

    call boundary__field(gkl,nx,bc)

    f3 = c*delt*gfac/delx
    do i=1,nx
       gkl(2,i) = gkl(2,i)+f3*(-gkl(6,i-1)+gkl(6,i))
       gkl(3,i) = gkl(3,i)-f3*(-gkl(5,i-1)+gkl(5,i))
    enddo

    !solve  < by & bz;  bx = const >
    call cgm(gf,gkl,nx,c,delx,delt,gfac,bc)
    call boundary__field(gf,nx,bc)

    !solve  < ex, ey & ez >
    do i=1,nx
       gf(4,i) = gkl(4,i)   ! if (Ex=0) gf(4,i)=0.
    enddo
    do i=1,nx+bc
       gf(5,i) = gkl(5,i)-f3*(-gf(3,i)+gf(3,i+1))
       gf(6,i) = gkl(6,i)+f3*(-gf(2,i)+gf(2,i+1))
    enddo

    call boundary__field(gf,nx,bc)

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


  subroutine ele_cur(uj,up,gp,np,nx,nsp,np2,bc,q,c,delx,delt)

    integer, intent(in)  :: np, nx, nsp, bc
    integer, intent(in)  :: np2(1:nx+bc,nsp)
    real(8), intent(in)  :: q(nsp), c, delx, delt
    real(8), intent(in)  :: up(4,np,1:nx+bc,nsp), gp(4,np,1:nx+bc,nsp)
    real(8), intent(out) :: uj(3,-1:nx+2)
    integer :: ii, i, isp, ih, i1 ,i2
    real(8) :: x2, xh, xr, qvx1, qvx2, idelt, idelx, gam
    real(8) :: dx, dxm

    !memory clear
    uj(1:3,-1:nx+2) = 0.0D0

    !------ Charge Conservation Method for Jx ---------!
    !----  Zigzag scheme (Umeda et al., CPC, 2003) ----!
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
             gam = 1./dsqrt(1.0+(+gp(2,ii,i,isp)*gp(2,ii,i,isp) &
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


  subroutine cgm(gb,gl,nx,c,delx,delt,gfac,bc)

    !-----------------------------------------------------------------------
    !  #  conjugate gradient method 
    !-----------------------------------------------------------------------

    integer, intent(in)    :: nx, bc
    real(8), intent(in)    :: c, delx, delt, gfac
    real(8), intent(in)    :: gl(6,0:nx+1)
    real(8), intent(inout) :: gb(6,0:nx+1)
    integer :: ite_max = 100 ! maximum number of interation
    integer :: i, l, ite
    real(8) :: err = 1d-6 
    real(8) :: f1, f2, eps, sumr, sum, sum1, sum2, av, bv
    real(8) :: x(0:nx+1), b(0:nx+1), r(0:nx+1), p(0:nx+1), ap(0:nx+1)


    do l=2,3

       ! initial guess
       ite = 0
       f1 = 2.0+(delx/(c*delt*gfac))**2
       f2 = (delx/(c*delt*gfac))**2
       sum = 0.0
       do i=1,nx
          x(i) = gb(l,i)
          b(i) = f2*gl(l,i)
          sum = sum+b(i)*b(i)
       enddo
       eps = dsqrt(sum)*err
       
       if(bc == 0)then
          !periodic condition
          x(0)    = x(nx)
          x(nx+1) = x(1)
       else if(bc == -1)then
          !reflective condition
          x(0)  = x(2)
          x(nx+1) = x(nx-1)
       else
          write(*,*)'choose bc=0 (periodic) or bc=-1 (reflective)'
       endif

       do i=1,nx
          r(i) = b(i)+x(i-1)-f1*x(i)+x(i+1)
          p(i) = r(i)
          sumr = sumr+r(i)*r(i)
       enddo

       if(dsqrt(sumr) > eps)then
       
          do while(sum > eps)
             
             ite = ite+1

             if(bc == 0)then
                !periodic condition
                p(0)    = p(nx)
                p(nx+1) = p(1)
             else if(bc == -1)then
                !reflective condition
                p(0)  = p(2)
                p(nx+1) = p(nx-1)
             else
                write(*,*)'choose bc=0 (periodic) or bc=-1 (reflective)'
                stop
             endif
       
             do i=1,nx
                ap(i) = -p(i-1)+f1*p(i)-p(i+1)
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
             
             sum = dsqrt(sumr)
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


!!$  subroutine ele_cur2(uj,up,np,nx,nsp,np2,bc,q,c,delx)
!!$
!!$    integer, intent(in)  :: np, nx, nsp, bc
!!$    integer, intent(in)  :: np2(1:nx+bc,nsp)
!!$    real(8), intent(in)  :: q(nsp), c, delx
!!$    real(8), intent(in)  :: up(4,np,1:nx+bc,nsp)
!!$    real(8), intent(out) :: uj(3,-1:nx+2)
!!$    integer :: ii, i, isp, ih
!!$    real(8) :: dx, dxm, gam, idelx
!!$
!!$    !memory clear
!!$    uj(1:3,-1:nx+2) = 0.0D0
!!$
!!$    !caluculate erectric current density
!!$    do isp=1,nsp
!!$       do i=1,nx+bc
!!$          do ii=1,np2(i,isp)
!!$             gam = 1./dsqrt(1.0+(+up(2,ii,i,isp)*up(2,ii,i,isp) &
!!$                                 +up(3,ii,i,isp)*up(3,ii,i,isp) &
!!$                                 +up(4,ii,i,isp)*up(4,ii,i,isp) &
!!$                                )/(c*c))
!!$
!!$             dx = up(1,ii,i,isp)-i
!!$             dxm = 1.-dx
!!$             uj(1,i  ) = uj(1,i  )+q(isp)*up(2,ii,i,isp)*gam*dxm
!!$             uj(1,i+1) = uj(1,i+1)+q(isp)*up(2,ii,i,isp)*gam*dx 
!!$
!!$             ih = floor(up(1,ii,i,isp)-0.5)
!!$             dx = up(1,ii,i,isp)-0.5-ih
!!$             dxm = 1.-dx
!!$
!!$             uj(2,ih  ) = uj(2,ih  )+q(isp)*up(3,ii,i,isp)*gam*dxm
!!$             uj(3,ih  ) = uj(3,ih  )+q(isp)*up(4,ii,i,isp)*gam*dxm
!!$             uj(2,ih+1) = uj(2,ih+1)+q(isp)*up(3,ii,i,isp)*gam*dx 
!!$             uj(3,ih+1) = uj(3,ih+1)+q(isp)*up(4,ii,i,isp)*gam*dx 
!!$          enddo
!!$       enddo
!!$    enddo
!!$
!!$    idelx = 1.D0/delx
!!$    if(bc == 0)then
!!$       do i=-1,nx+2
!!$          uj(1:3,i) = uj(1:3,i)*idelx
!!$       enddo
!!$    else if(bc == -1)then
!!$       i=1
!!$       uj(1,i) = uj(1,i)*2.*idelx
!!$       do i=2,nx-1
!!$          uj(1,i) = uj(1,i)*idelx
!!$       enddo
!!$       i=nx
!!$       uj(1,i) = uj(1,i)*2.*idelx
!!$
!!$       do i=0,nx
!!$          uj(2,i) = uj(2,i)*idelx
!!$          uj(3,i) = uj(3,i)*idelx
!!$       enddo
!!$    else
!!$       write(*,*)'choose bc=0 (periodic) or bc=-1 (reflective)'
!!$       stop
!!$    endif
!!$
!!$  end subroutine ele_cur2


end module field
