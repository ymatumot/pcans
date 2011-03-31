module field

  implicit none

  private

  public :: field__fdtd_i


contains

  
  subroutine field__fdtd_i(uf,up,gp,                                      &
                           np,nsp,np2,nxgs,nxge,nxs,nxe,nxs1,nxe1,bc,bcp, &
                           q,c,delx,delt,gfac,                            &
                           nup,ndown,nroot,nproc,nrank,mnpr,opsum,nstat,ncomw,nerr)

    use boundary, only : boundary__field, boundary__curre,  boundary__particle
 
    integer, intent(in)    :: np, nsp, nxgs, nxge, nxs, nxe, nxs1, nxe1, bc, bcp
    integer, intent(in)    :: np2(nxs:nxe+bcp,nsp) 
    integer, intent(in)    :: nup, ndown, nroot, nproc, nrank, opsum, mnpr, ncomw
    integer, intent(inout) :: nerr, nstat(:)
    real(8), intent(in)    :: q(nsp), c, delx, delt, gfac
    real(8), intent(in)    :: gp(4,np,nxs:nxe+bcp,nsp)
    real(8), intent(inout) :: up(4,np,nxs:nxe+bcp,nsp)
    real(8), intent(inout) :: uf(6,nxs1:nxe1)
    integer                    :: ii, i, isp
    integer, save              :: flag
    real(8)                    :: uj(3,nxs1-1:nxe1+1), gkl(6,nxs1:nxe1)
    real(8), save, allocatable :: gf(:,:)
    real(8)                    :: pi, f1, f2, f3, rotb2, rotb3, rote2, rote3

    pi = 4.0*atan(1.0)

    if(flag /=1)then
       allocate(gf(6,nxs1:nxe1))
       gf(1:3,nxs1:nxe1) = 0.0
       flag=1
    endif

    !position at n+1/2
!!$    do isp=1,nsp
!!$       do i=nxs,nxe+bcp
!!$          do ii=1,np2(i,isp)
!!$             up(1,ii,i,isp) = 0.5*(up(1,ii,i,isp)+gp(1,ii,i,isp))
!!$             up(2:4,ii,i,isp) = gp(2:4,ii,i,isp)
!!$          enddo
!!$       enddo
!!$    enddo
!!$    call ele_cur2(uj,up, &
!!$                  np,nsp,np2,nxs,nxe,nxs1,nxe1,bcp,q,c)
    call ele_cur(uj,up,gp,                                      &
                 np,nsp,np2,nxgs,nxge,nxs,nxe,nxs1,nxe1,bc,bcp, &
                 q,c,delx,delt,                                 &
                 nroot,nproc,nrank)
    call boundary__curre(uj,nxs,nxe,nxs1,nxe1,bc, &
                         nup,ndown,nroot,nproc,nrank,mnpr,nstat,ncomw,nerr)

    !calculation
    !< gk =  (c*delt)*rot(b) - (4*pi*delt)*j >
    !< gl =  (c*delt)*rot(e) >
    f1 = c*delt/delx
    f2 = 4.0*pi*delt
    !By,z
    do i=nxs,nxe
       rote3 = -uf(5,i-1)+uf(5,i)
       rote2 = +uf(6,i-1)-uf(6,i)
       gkl(2,i) = -f1*rote2
       gkl(3,i) = -f1*rote3
    enddo
    !Ex
    do i=nxs,nxe
       gkl(4,i) = -f2*uj(1,i)
    enddo
    !Ey,z
    do i=nxs,nxe+bcp
       rotb3 = -uf(2,i)+uf(2,i+1)
       rotb2 = +uf(3,i)-uf(3,i+1)
       gkl(5,i) = f1*rotb2-f2*uj(2,i)
       gkl(6,i) = f1*rotb3-f2*uj(3,i)
    enddo

    call boundary__field(gkl,                      &
                         nxs,nxe,nxs1,nxe1,bc,bcp, &
                         nup,ndown,nroot,nproc,nrank,mnpr,nstat,ncomw,nerr)

    f3 = c*delt*gfac/delx
    do i=nxs,nxe
       gkl(2,i) = gkl(2,i)+f3*(-gkl(6,i-1)+gkl(6,i))
       gkl(3,i) = gkl(3,i)-f3*(-gkl(5,i-1)+gkl(5,i))
    enddo

    !solve  < by & bz;  bx = const >
    call cgm(gf,gkl,              &
             nxs,nxe,nxs1,nxe1,   &
             c,delx,delt,gfac,bc, &
             nup,ndown,nroot,nproc,nrank,mnpr,opsum,nstat,ncomw,nerr)
    
    call boundary__field(gf,                       &
                         nxs,nxe,nxs1,nxe1,bc,bcp, &
                         nup,ndown,nroot,nproc,nrank,mnpr,nstat,ncomw,nerr)

    !solve  < ex, ey & ez >
    do i=nxs,nxe
       gf(4,i) = gkl(4,i)
    enddo
    do i=nxs,nxe+bcp
       gf(5,i) = gkl(5,i)-f3*(-gf(3,i)+gf(3,i+1))
       gf(6,i) = gkl(6,i)+f3*(-gf(2,i)+gf(2,i+1))
    enddo

    call boundary__field(gf,                       &
                         nxs,nxe,nxs1,nxe1,bc,bcp, &
                         nup,ndown,nroot,nproc,nrank,mnpr,nstat,ncomw,nerr)

    !Update
    uf(1:6,nxs1:nxe1) = uf(1:6,nxs1:nxe1)+gf(1:6,nxs1:nxe1)

    do isp=1,nsp
       do i=nxs,nxe+bcp
          do ii=1,np2(i,isp)
             up(1,ii,i,isp) = gp(1,ii,i,isp)
             up(2,ii,i,isp) = gp(2,ii,i,isp)
             up(3,ii,i,isp) = gp(3,ii,i,isp)
             up(4,ii,i,isp) = gp(4,ii,i,isp)
          enddo
       enddo
    enddo

  end subroutine field__fdtd_i


  subroutine ele_cur(uj,up,gp,                                      &
                     np,nsp,np2,nxgs,nxge,nxs,nxe,nxs1,nxe1,bc,bcp, &
                     q,c,delx,delt,                                 &
                     nroot,nproc,nrank)

    integer, intent(in)  :: np, nsp, nxgs, nxge, nxs, nxe, nxs1, nxe1, bc, bcp
    integer, intent(in)  :: np2(nxs:nxe+bcp,nsp)
    integer, intent(in)  :: nroot, nproc, nrank
    real(8), intent(in)  :: q(nsp), c, delx, delt
    real(8), intent(in)  :: up(4,np,nxs:nxe+bcp,nsp), gp(4,np,nxs:nxe+bcp,nsp)
    real(8), intent(out) :: uj(3,nxs1-1:nxe1+1)
    integer :: ii, i, isp, ih, i1 ,i2
    real(8) :: x2, xh, xr, qvx1, qvx2, idelt, idelx, gam
    real(8) :: dx, dxm

    !memory clear
    uj(1:3,nxs1-1:nxe1+1) = 0.0D0

    !------ Charge Conservation Method for Jx ---------!
    !----  Zigzag scheme (Umeda et al., CPC, 2003) ----!
    idelt = 1.D0/delt
    idelx = 1.D0/delx
    do isp=1,nsp
       do i=nxs,nxe+bcp
          do ii=1,np2(i,isp)

             x2 = gp(1,ii,i,isp)

             !reflective boundary condition in x
             if(bc == -1)then
                if(x2 < nxgs)then
                   x2  = 2.*nxgs-x2
                endif
                if(x2 > nxge)then
                   x2  = 2.*nxge-x2
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
             uj(1,i1+1) = uj(1,i1+1)+qvx1*idelx
             uj(1,i2+1) = uj(1,i2+1)+qvx2*idelx

             !Jy and Jz
             gam = 1./dsqrt(1.0+(+gp(2,ii,i,isp)*gp(2,ii,i,isp) &
                                 +gp(3,ii,i,isp)*gp(3,ii,i,isp) &
                                 +gp(4,ii,i,isp)*gp(4,ii,i,isp) &
                                )/(c*c))

             ih = floor(xh-0.5)
             dx = xh-0.5-ih
             dxm = 1.-dx

             uj(2,ih  ) = uj(2,ih  )+q(isp)*gp(3,ii,i,isp)*gam*dxm*idelx
             uj(3,ih  ) = uj(3,ih  )+q(isp)*gp(4,ii,i,isp)*gam*dxm*idelx
             uj(2,ih+1) = uj(2,ih+1)+q(isp)*gp(3,ii,i,isp)*gam*dx *idelx
             uj(3,ih+1) = uj(3,ih+1)+q(isp)*gp(4,ii,i,isp)*gam*dx *idelx
          enddo
       enddo
    enddo

    if(bc == -1)then
       if(nrank == nroot)then
          i=nxs
          uj(1,i) = uj(1,i)*2.*idelx
       endif
       if(nrank == nproc-1)then
          i=nxe
          uj(1,i) = uj(1,i)*2.*idelx
       endif
    endif

  end subroutine ele_cur


  subroutine cgm(gb,gl,               &
                 nxs,nxe,nxs1,nxe1,   &
                 c,delx,delt,gfac,bc, &
                 nup,ndown,nroot,nproc,nrank,mnpr,opsum,nstat,ncomw,nerr)

    !-----------------------------------------------------------------------
    !  #  conjugate gradient method 
    !  #  this routine will be stoped after itaration number = ite_max
    !-----------------------------------------------------------------------

    integer, intent(in)    :: nxs, nxe, nxs1, nxe1, bc
    integer, intent(in)    :: nup, ndown, nroot, nproc, nrank, mnpr, opsum, ncomw
    integer, intent(inout) :: nerr, nstat(:)
    real(8), intent(in)    :: c, delx, delt, gfac
    real(8), intent(in)    :: gl(6,nxs1:nxe1)
    real(8), intent(inout) :: gb(6,nxs1:nxe1)
    integer, parameter :: ite_max = 100 ! maximum number of interation
    integer            :: i, l, ite
    real(8), parameter :: err = 1d-6 
    real(8)            :: f1, f2, eps, sumr, sum, sum1, sum2, av, bv
    real(8)            :: sumr_g, sum_g, sum1_g, sum2_g, bff_snd(2), bff_rcv(2)
    real(8)            :: x(nxs1:nxe1), b(nxs1:nxe1), r(nxs1:nxe1), p(nxs1:nxe1), ap(nxs1:nxe1)

    do l=2,3

       ! initial guess
       ite = 0
       f1 = 2.0+(delx/(c*delt*gfac))**2
       f2 = (delx/(c*delt*gfac))**2
       sum = 0.0
       do i=nxs,nxe
          x(i) = gb(l,i)
          b(i) = f2*gl(l,i)
          sum = sum+b(i)*b(i)
       enddo

       call MPI_ALLREDUCE(sum,sum_g,1,mnpr,opsum,ncomw,nerr)

       eps = dsqrt(sum_g)*err

       !boundary condition of x
       call MPI_SENDRECV(x(nxe),1,mnpr,nup,100,    &
                         x(nxs1),1,mnpr,ndown,100, &
                         ncomw,nstat,nerr)
       call MPI_SENDRECV(x(nxs),1,mnpr,ndown,101,  &
                         x(nxe1),1,mnpr,nup,101,   &
                         ncomw,nstat,nerr)
       if(bc == -1)then !reflective condition
          if(nrank == nroot)then
             x(nxs1) = +x(nxs+1)
          endif
          if(nrank == nproc-1)then
             x(nxe1) = +x(nxe-1)
          endif
       endif
          
       do i=nxs,nxe
          r(i) = b(i)+x(i-1)-f1*x(i)+x(i+1)
          p(i) = r(i)
          sumr = sumr+r(i)*r(i)
       enddo
       call MPI_ALLREDUCE(sumr,sumr_g,1,mnpr,opsum,ncomw,nerr)

       if(dsqrt(sumr_g) > eps)then
       
          do while(sum_g > eps)
             
             ite = ite+1

             !boundary condition of p
             call MPI_SENDRECV(p(nxe),1,mnpr,nup,100,    &
                               p(nxs1),1,mnpr,ndown,100, &
                               ncomw,nstat,nerr)
             call MPI_SENDRECV(p(nxs),1,mnpr,ndown,101,  &
                               p(nxe1),1,mnpr,nup,101,   &
                               ncomw,nstat,nerr)
             if(bc == -1)then !reflective condition
                if(nrank == nroot)then
                   p(nxs1) = +p(nxs+1)
                endif
                if(nrank == nproc-1)then
                   p(nxe1) = +p(nxe-1)
                endif
             endif

             do i=nxs,nxe
                ap(i) = -p(i-1)+f1*p(i)-p(i+1)
             enddo
             sumr = 0.0
             sum2 = 0.0
             do i=nxs,nxe
                sumr = sumr+r(i)*r(i)
                sum2 = sum2+p(i)*ap(i)
             enddo

             bff_snd(1) = sumr
             bff_snd(2) = sum2
             call MPI_ALLREDUCE(bff_snd,bff_rcv,2,mnpr,opsum,ncomw,nerr)
             sumr_g = bff_rcv(1)
             sum2_g = bff_rcv(2)

             av = sumr_g/sum2_g
             
             x(nxs:nxe) = x(nxs:nxe)+av*p(nxs:nxe)
             r(nxs:nxe) = r(nxs:nxe)-av*ap(nxs:nxe)
             
             sum_g = dsqrt(sumr_g)
             if(ite >= ite_max) then
                write(6,*)'********** stop at cgm after ite_max **********'
                stop
             endif
             
             sum1 = 0.0
             do i=nxs,nxe
                sum1 = sum1+r(i)*r(i)
             end do
             call MPI_ALLREDUCE(sum1,sum1_g,1,mnpr,opsum,ncomw,nerr)
             bv = sum1_g/sumr_g
             
             p(nxs:nxe) = r(nxs:nxe)+bv*p(nxs:nxe)
             
          enddo
       endif

       gb(l,nxs:nxe) = x(nxs:nxe)

    end do
    
  end subroutine cgm


  subroutine ele_cur2(uj,up, &
                      np,nsp,np2,nxs,nxe,nxs1,nxe1,bcp,q,c)

    integer, intent(in)  :: np, nsp, nxs, nxe, nxs1, nxe1, bcp
    integer, intent(in)  :: np2(nxs:nxe+bcp,nsp)
    real(8), intent(in)  :: q(nsp), c
    real(8), intent(in)  :: up(4,np,nxs:nxe+bcp,nsp)
    real(8), intent(out) :: uj(3,nxs1-1:nxe1+1)
    integer :: ii, i, isp, ih
    real(8) :: dx, dxm, dxh, dxmh, gam

    !memory clear
    uj(1:3,nxs1-1:nxe1+1) = 0.0D0

    !caluculate erectric current density
    do isp=1,nsp
       do i=nxs,nxe+bcp
          do ii=1,np2(i,isp)
             gam = 1./dsqrt(1.0+(+up(2,ii,i,isp)*up(2,ii,i,isp) &
                                 +up(3,ii,i,isp)*up(3,ii,i,isp) &
                                 +up(4,ii,i,isp)*up(4,ii,i,isp) &
                                )/(c*c))

             dx   = up(1,ii,i,isp)-i
             dxm  = 1.-dx

             ih = floor(up(1,ii,i,isp)+0.5)
             dxh  = up(1,ii,i,isp)+0.5-ih
             dxmh = 1.-dxh

             uj(1   ,i) = uj(1   ,i)+q(isp)*up(2,ii,i,isp)*gam*dxm
             uj(1 ,i+1) = uj(1 ,i+1)+q(isp)*up(2,ii,i,isp)*gam*dx 

             uj(2,ih-1) = uj(2,ih-1)+q(isp)*up(3,ii,i,isp)*gam*dxmh
             uj(3,ih-1) = uj(3,ih-1)+q(isp)*up(4,ii,i,isp)*gam*dxmh

             uj(2  ,ih) = uj(2  ,ih)+q(isp)*up(3,ii,i,isp)*gam*dxh
             uj(3  ,ih) = uj(3  ,ih)+q(isp)*up(4,ii,i,isp)*gam*dxh
          enddo
       enddo
    enddo

  end subroutine ele_cur2


end module field
