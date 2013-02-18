module field

  implicit none

  private

  public :: field__fdtd_i


contains

  
  subroutine field__fdtd_i(uf,up,gp,                           &
                           np,nsp,np2,nxs,nxe,nys,nye,nsfo,bc, &
                           q,c,delx,delt,gfac,                 &
                           nup,ndown,mnpr,opsum,nstat,ncomw,nerr)

    use boundary, only : boundary__field, boundary__curre,  boundary__particle
 
    integer, intent(in)    :: np, nsp, nxs, nxe, nys, nye, nsfo, bc
    integer, intent(in)    :: np2(nys:nye,nsp)
    integer, intent(in)    :: nup, ndown, opsum, mnpr, ncomw
    integer, intent(inout) :: nerr, nstat(:)
    real(8), intent(in)    :: q(nsp), c, delx, delt, gfac
    real(8), intent(in)    :: gp(5,np,nys:nye,nsp)
    real(8), intent(inout) :: up(5,np,nys:nye,nsp)
    real(8), intent(inout) :: uf(6,nxs-2:nxe+2,nys-2:nye+2)
    logical, save              :: lflag=.true.
    integer                    :: ii, i, j, isp
    real(8)                    :: pi, f1, f2, f3
    real(8)                    :: uj(3,nxs-2:nxe+2,nys-2:nye+2)
    real(8)                    :: gkl(6,nxs-2:nxe+2,nys-2:nye+2)
    real(8), save, allocatable :: gf(:,:,:)

    pi = 4.0*atan(1.0)

    if(lflag)then
       allocate(gf(6,nxs-2:nxe+2,nys-2:nye+2))
       gf(1:6,nxs-2:nxe+2,nys-2:nye+2) = 0.0
       lflag = .false.
    endif

    call ele_cur(uj,up,gp, &
                 np,nsp,np2,nxs,nxe,nys,nye,bc,nsfo,q,c,delx,delt)
    call boundary__curre(uj,nxs,nxe,nys,nye,bc, &
                         nup,ndown,mnpr,nstat,ncomw,nerr)

    !calculation
    !< gkl(1:3) =  (c*delt)*rot(e) >
    !< gkl(4:6) =  (c*delt)*rot(b) - (4*pi*delt)*j >
    f1 = c*delt/delx
    f2 = 4.0*pi*delt
    do j=nys,nye
    do i=nxs,nxe+bc
       gkl(1,i,j) = -f1*(+(-uf(6,i,j-1)+uf(6,i,j)))
       gkl(2,i,j) = -f1*(-(-uf(6,i-1,j)+uf(6,i,j)))
       gkl(3,i,j) = -f1*(-(-uf(4,i,j-1)+uf(4,i,j))+(-uf(5,i-1,j)+uf(5,i,j)))
       gkl(4,i,j) = +f1*(+(-uf(3,i,j)+uf(3,i,j+1)))-f2*uj(1,i,j)
       gkl(5,i,j) = +f1*(-(-uf(3,i,j)+uf(3,i+1,j)))-f2*uj(2,i,j)
       gkl(6,i,j) = +f1*(-(-uf(1,i,j)+uf(1,i,j+1))+(-uf(2,i,j)+uf(2,i+1,j)))-f2*uj(3,i,j)
    enddo
    enddo
    if(bc == -1)then
       i=nxe
       do j=nys,nye
          gkl(2,i,j) = -f1*(-(-uf(6,i-1,j)+uf(6,i,j)))
          gkl(3,i,j) = -f1*(-(-uf(4,i,j-1)+uf(4,i,j))+(-uf(5,i-1,j)+uf(5,i,j)))
          gkl(4,i,j) = +f1*(+(-uf(3,i,j)+uf(3,i,j+1)))-f2*uj(1,i,j)
       enddo
    endif

    call boundary__field(gkl,                &
                         nxs,nxe,nys,nye,bc, &
                         nup,ndown,mnpr,nstat,ncomw,nerr)

    f3 = c*delt*gfac/delx
    do j=nys,nye
    do i=nxs,nxe+bc
       gkl(1,i,j) = gkl(1,i,j)-f3*(-gkl(6,i,j-1)+gkl(6,i,j))
       gkl(2,i,j) = gkl(2,i,j)+f3*(-gkl(6,i-1,j)+gkl(6,i,j))
       gkl(3,i,j) = gkl(3,i,j)-f3*(+gkl(4,i,j-1)-gkl(4,i,j) &
                                   -gkl(5,i-1,j)+gkl(5,i,j))
    enddo
    enddo
    if(bc == -1)then
       i=nxe
       do j=nys,nye
          gkl(2,i,j) = gkl(2,i,j)+f3*(-gkl(6,i-1,j)+gkl(6,i,j))
          gkl(3,i,j) = gkl(3,i,j)-f3*(+gkl(4,i,j-1)-gkl(4,i,j) &
                                      -gkl(5,i-1,j)+gkl(5,i,j))
       enddo
    endif

    !solve  < bx, by & bz >
    call cgm(gf,gkl,             &
             nxs,nxe,nys,nye,bc, &
             c,delx,delt,gfac,   &
             nup,ndown,mnpr,opsum,nstat,ncomw,nerr)
    
    call boundary__field(gf,                 &
                         nxs,nxe,nys,nye,bc, &
                         nup,ndown,mnpr,nstat,ncomw,nerr)

    !solve  < ex, ey & ez >
    do j=nys,nye
    do i=nxs,nxe+bc
       gf(4,i,j) = gkl(4,i,j)+f3*(-gf(3,i,j)+gf(3,i,j+1))
       gf(5,i,j) = gkl(5,i,j)-f3*(-gf(3,i,j)+gf(3,i+1,j))
       gf(6,i,j) = gkl(6,i,j)+f3*(-gf(2,i,j)+gf(2,i+1,j) &
                                  +gf(1,i,j)-gf(1,i,j+1))
    enddo
    enddo
    if(bc == -1)then
       i=nxe
       do j=nys,nye
          gf(4,i,j) = gkl(4,i,j)+f3*(-gf(3,i,j)+gf(3,i,j+1))
       enddo
    endif

    call boundary__field(gf,                 &
                         nxs,nxe,nys,nye,bc, &
                         nup,ndown,mnpr,nstat,ncomw,nerr)

    !===== Update fields and particles ======
    uf(1:6,nxs-2:nxe+2,nys-2:nye+2) = uf(1:6,nxs-2:nxe+2,nys-2:nye+2) &
                                     +gf(1:6,nxs-2:nxe+2,nys-2:nye+2)

    do isp=1,nsp
       do j=nys,nye
          do ii=1,np2(j,isp)
             up(1,ii,j,isp) = gp(1,ii,j,isp)
             up(2,ii,j,isp) = gp(2,ii,j,isp)
             up(3,ii,j,isp) = gp(3,ii,j,isp)
             up(4,ii,j,isp) = gp(4,ii,j,isp)
             up(5,ii,j,isp) = gp(5,ii,j,isp)
          enddo
       enddo
    enddo

  end subroutine field__fdtd_i


  subroutine ele_cur(uj,up,gp, &
                     np,nsp,np2,nxs,nxe,nys,nye,bc,nsfo,q,c,delx,delt)

    use shape_function, only : sf

    integer, intent(in)  :: np, nsp, nxs, nxe, nys, nye, bc, nsfo
    integer, intent(in)  :: np2(nys:nye,nsp)
    real(8), intent(in)  :: q(nsp), c, delx, delt
    real(8), intent(in)  :: up(5,np,nys:nye,nsp), gp(5,np,nys:nye,nsp)
    real(8), intent(out) :: uj(3,nxs-2:nxe+2,nys-2:nye+2)
    
    integer            :: ii, j, isp, i1, j1, i2, j2, iinc, jinc, ip, jp
    real(8), parameter :: fac = 1.D0/3.D0
    real(8)            :: idelx, idelt, gamp
    real(8)            :: s0(-2:2,2), s1(-2:2,2), ds(-2:2,2), ujp(3,-2:3,-2:3)

    uj(1:3,nxs-2:nxe+2,nys-2:nye+2) = 0.D0

    idelt = 1.D0/delt
    idelx = 1.D0/delx

    !--------------Charge Conservation Method -------------!
    !---- Density Decomposition (Esirkepov, CPC, 2001) ----!
    do isp=1,nsp
       do j=nys,nye
          do ii=1,np2(j,isp)

             ujp(1:3,-2:3,-2:3) = 0.D0

             i1 = int(up(1,ii,j,isp)*idelx)
             i2 = int(gp(1,ii,j,isp)*idelx)
             j2 = int(gp(2,ii,j,isp)*idelx)
             iinc = i2-i1
             jinc = j2-j

             gamp = 1./dsqrt(1.+(+gp(3,ii,j,isp)*gp(3,ii,j,isp) &
                                 +gp(4,ii,j,isp)*gp(4,ii,j,isp) &
                                 +gp(5,ii,j,isp)*gp(5,ii,j,isp))/(c*c) )

             s0(-2:2,1) = sf(i1,up(1,ii,j,isp)*idelx-0.5,nsfo)
             s0(-2:2,2) = sf(j ,up(2,ii,j,isp)*idelx-0.5,nsfo)

             s1(-2:2,1) = sf(i2,gp(1,ii,j,isp)*idelx-0.5,nsfo)
             s1(-2:2,2) = sf(j2,gp(2,ii,j,isp)*idelx-0.5,nsfo)

             ds(-2:2,1) = cshift(s1(-2:2,1),-iinc,1)-s0(-2:2,1)
             ds(-2:2,2) = cshift(s1(-2:2,2),-jinc,1)-s0(-2:2,2)

             do jp=-1+min(jinc,0),1+max(jinc,0)
             do ip=-1+min(iinc,0),1+max(iinc,0)
                ujp(1,ip+1,jp) = ujp(1,ip,jp) &
                                -q(isp)*delx*idelt*ds(ip,1)*(s0(jp,2)+0.5*ds(jp,2))

                ujp(2,ip,jp+1) = ujp(2,ip,jp) &
                                -q(isp)*delx*idelt*ds(jp,2)*(s0(ip,1)+0.5*ds(ip,1))

                ujp(3,ip,jp) = +q(isp)*gp(5,ii,j,isp)*gamp                &
                               *(+s0(ip,1)*s0(jp,2)+0.5*ds(ip,1)*s0(jp,2) &
                                 +0.5*s0(ip,1)*ds(jp,2)+fac*ds(ip,1)*ds(jp,2))
             enddo
             enddo

             uj(1:3,i1-2:i1+2,j-2:j+2) = uj(1:3,i1-2:i1+2,j-2:j+2)+ujp(1:3,-2:2,-2:2)

          enddo
       enddo
    enddo

  end subroutine ele_cur


  subroutine cgm(gb,gkl,             &
                 nxs,nxe,nys,nye,bc, &
                 c,delx,delt,gfac,   &
                 nup,ndown,mnpr,opsum,nstat,ncomw,nerr)

    !-----------------------------------------------------------------------
    !  #  conjugate gradient method 
    !  #  this routine will be stopped after iteration number = ite_max
    !-----------------------------------------------------------------------

    integer, intent(in)    :: nxs, nxe, nys, nye, bc
    integer, intent(in)    :: nup, ndown, mnpr, opsum, ncomw
    integer, intent(inout) :: nerr, nstat(:)
    real(8), intent(in)    :: c, delx, delt, gfac
    real(8), intent(in)    :: gkl(6,nxs-2:nxe+2,nys-2:nye+2)
    real(8), intent(inout) :: gb(6,nxs-2:nxe+2,nys-2:nye+2)
    integer, parameter :: ite_max = 100 ! maximum number of iteration
    integer            :: i, ii, j, l, ite
    real(8), parameter :: err = 1d-6 
    real(8)            :: f1, f2, eps, sumr, sum, sum1, sum2, av, bv
    real(8)            :: sumr_g, sum_g, sum1_g, sum2_g
    real(8)            :: x(nxs-1:nxe+1,nys-1:nye+1), b(nxs:nxe,nys:nye)
    real(8)            :: r(nxs:nxe,nys:nye), p(nxs-1:nxe+1,nys-1:nye+1)
    real(8)            :: ap(nxs:nxe,nys:nye)
    real(8)            :: bff_snd(nxe-nxs+1), bff_rcv(nxe-nxs+1)

    do l=1,1

       ! initial guess
       ite = 0
       f2 = (delx/(c*delt*gfac))**2
       sum = 0.0
       do j=nys,nye
       do i=nxs,nxe+bc
          x(i,j) = gb(l,i,j)
          b(i,j) = f2*gkl(l,i,j)
          sum = sum+b(i,j)*b(i,j)
       enddo
       enddo

       call MPI_ALLREDUCE(sum,sum_g,1,mnpr,opsum,ncomw,nerr)

       eps = dsqrt(sum_g)*err

       !------ boundary condition of x ------
       do i=nxs,nxe+bc
          ii = i-nxs+1
          bff_snd(ii) = x(i,nys)
       enddo

       call MPI_SENDRECV(bff_snd(1),nxe+bc-nxs+1,mnpr,ndown,101, &
                         bff_rcv(1),nxe+bc-nxs+1,mnpr,nup  ,101, &
                         ncomw,nstat,nerr)

       do i=nxs,nxe+bc
          ii = i-nxs+1
          x(i,nye+1) = bff_rcv(ii)
       enddo

       do i=nxs,nxe+bc
          ii = i-nxs+1
          bff_snd(ii) = x(i,nye)
       enddo

       call MPI_SENDRECV(bff_snd(1),nxe+bc-nxs+1,mnpr,nup  ,100, &
                         bff_rcv(1),nxe+bc-nxs+1,mnpr,ndown,100, &
                         ncomw,nstat,nerr)

       do i=nxs,nxe+bc
          ii = i-nxs+1
          x(i,nys-1) = bff_rcv(ii)
       enddo

       if(bc == 0)then
          x(nxs-1,nys-1:nye+1) = x(nxe,nys-1:nye+1)
          x(nxe+1,nys-1:nye+1) = x(nxs,nys-1:nye+1)
       else if(bc == -1)then
          x(nxs-1,nys-1:nye+1) = -x(nxs  ,nys-1:nye+1)
          x(nxe  ,nys-1:nye+1) = -x(nxe-1,nys-1:nye+1)
       else
          write(*,*)'choose bc=0 (periodic) or bc=-1 (reflective)'
          stop
       endif
       !------ end of -----

       f1 = 4.0+(delx/(c*delt*gfac))**2
       sumr = 0.0
       do j=nys,nye
       do i=nxs,nxe+bc
          r(i,j) = b(i,j)+x(i,j-1)                    &
                         +x(i-1,j)-f1*x(i,j)+x(i+1,j) &
                         +x(i,j+1)
          p(i,j) = r(i,j)
          sumr = sumr+r(i,j)*r(i,j)
       enddo
       enddo

       call MPI_ALLREDUCE(sumr,sumr_g,1,mnpr,opsum,ncomw,nerr)

       if(dsqrt(sumr_g) > eps)then
       
          do while(sum_g > eps)
             
             ite = ite+1

             !------boundary condition of p------
             do i=nxs,nxe+bc
                ii = i-nxs+1
                bff_snd(ii) = p(i,nys)
             enddo

             call MPI_SENDRECV(bff_snd(1),nxe+bc-nxs+1,mnpr,ndown,101, &
                               bff_rcv(1),nxe+bc-nxs+1,mnpr,nup  ,101, &
                               ncomw,nstat,nerr)
             do i=nxs,nxe+bc
                ii = i-nxs+1
                p(i,nye+1) = bff_rcv(ii)
             enddo
             
             do i=nxs,nxe+bc
                ii = i-nxs+1
                bff_snd(ii) = p(i,nye)
             enddo
             call MPI_SENDRECV(bff_snd(1),nxe+bc-nxs+1,mnpr,nup  ,100, &
                               bff_rcv(1),nxe+bc-nxs+1,mnpr,ndown,100, &
                               ncomw,nstat,nerr)
             do i=nxs,nxe+bc
                ii = i-nxs+1
                p(i,nys-1) = bff_rcv(ii)
             enddo

             if(bc == 0)then
                p(nxs-1,nys-1:nye+1) = p(nxe,nys-1:nye+1)
                p(nxe+1,nys-1:nye+1) = p(nxs,nys-1:nye+1)
             else if(bc == -1)then
                p(nxs-1,nys-1:nye+1) = -p(nxs  ,nys-1:nye+1)
                p(nxe  ,nys-1:nye+1) = -p(nxe-1,nys-1:nye+1)
             else
                write(*,*)'choose bc=0 (periodic) or bc=-1 (reflective)'
                stop
             endif
             !----- end of -----
       
             do j=nys,nye
             do i=nxs,nxe+bc
                ap(i,j) = -p(i,j-1)                    &
                          -p(i-1,j)+f1*p(i,j)-p(i+1,j) &
                          -p(i,j+1)
             enddo
             enddo

             sumr = 0.0
             sum2 = 0.0
             do j=nys,nye
             do i=nxs,nxe+bc
                sumr = sumr+r(i,j)*r(i,j)
                sum2 = sum2+p(i,j)*ap(i,j)
             enddo
             enddo

             bff_snd(1) = sumr
             bff_snd(2) = sum2
             call MPI_ALLREDUCE(bff_snd,bff_rcv,2,mnpr,opsum,ncomw,nerr)
             sumr_g = bff_rcv(1)
             sum2_g = bff_rcv(2)

             av = sumr_g/sum2_g
             
             x(nxs:nxe+bc,nys:nye) = x(nxs:nxe+bc,nys:nye)+av* p(nxs:nxe+bc,nys:nye)
             r(nxs:nxe+bc,nys:nye) = r(nxs:nxe+bc,nys:nye)-av*ap(nxs:nxe+bc,nys:nye)
             
             sum_g = dsqrt(sumr_g)
             if(ite >= ite_max) then
                write(6,*)'********** stop at cgm after ite_max **********'
                stop
             endif
             
             sum1 = 0.0
             do j=nys,nye
             do i=nxs,nxe+bc
                sum1 = sum1+r(i,j)*r(i,j)
             enddo
             enddo
             call MPI_ALLREDUCE(sum1,sum1_g,1,mnpr,opsum,ncomw,nerr)
             bv = sum1_g/sumr_g
             
             p(nxs:nxe+bc,nys:nye) = r(nxs:nxe+bc,nys:nye)+bv*p(nxs:nxe+bc,nys:nye)
             
          enddo
       endif

       gb(l,nxs:nxe+bc,nys:nye) = x(nxs:nxe+bc,nys:nye)

    end do

    do l=2,3

       ! initial guess
       ite = 0
       f2 = (delx/(c*delt*gfac))**2
       sum = 0.0
       do j=nys,nye
       do i=nxs,nxe
          x(i,j) = gb(l,i,j)
          b(i,j) = f2*gkl(l,i,j)
          sum = sum+b(i,j)*b(i,j)
       enddo
       enddo

       call MPI_ALLREDUCE(sum,sum_g,1,mnpr,opsum,ncomw,nerr)

       eps = dsqrt(sum_g)*err

       !------ boundary condition of x ------
       do i=nxs,nxe
          ii = i-nxs+1
          bff_snd(ii) = x(i,nys)
       enddo
       call MPI_SENDRECV(bff_snd(1),nxe-nxs+1,mnpr,ndown,101, &
                         bff_rcv(1),nxe-nxs+1,mnpr,nup  ,101, &
                         ncomw,nstat,nerr)
       do i=nxs,nxe
          ii = i-nxs+1
          x(i,nye+1) = bff_rcv(ii)
       enddo

       do i=nxs,nxe
          ii = i-nxs+1
          bff_snd(ii) = x(i,nye)
       enddo
       call MPI_SENDRECV(bff_snd(1),nxe-nxs+1,mnpr,nup  ,100, &
                         bff_rcv(1),nxe-nxs+1,mnpr,ndown,100, &
                         ncomw,nstat,nerr)
       do i=nxs,nxe
          ii = i-nxs+1
          x(i,nys-1) = bff_rcv(ii)
       enddo

       if(bc == 0)then
          x(nxs-1,nys-1:nye+1) = x(nxe,nys-1:nye+1)
          x(nxe+1,nys-1:nye+1) = x(nxs,nys-1:nye+1)
       else if(bc == -1)then
          x(nxs-1,nys-1:nye+1) = x(nxs+1,nys-1:nye+1)
          x(nxe+1,nys-1:nye+1) = x(nxe-1,nys-1:nye+1)
       else
          write(*,*)'choose bc=0 (periodic) or bc=-1 (reflective)'
          stop
       endif
       !------ end of -----

       f1 = 4.0+(delx/(c*delt*gfac))**2
       do j=nys,nye
       do i=nxs,nxe
          r(i,j) = b(i,j)+x(i,j-1)                    &
                         +x(i-1,j)-f1*x(i,j)+x(i+1,j) &
                         +x(i,j+1)
          p(i,j) = r(i,j)
          sumr = sumr+r(i,j)*r(i,j)
       enddo
       enddo

       call MPI_ALLREDUCE(sumr,sumr_g,1,mnpr,opsum,ncomw,nerr)

       if(dsqrt(sumr_g) > eps)then
       
          do while(sum_g > eps)
             
             ite = ite+1

             !------boundary condition of p------
             do i=nxs,nxe
                ii = i-nxs+1
                bff_snd(ii) = p(i,nys)
             enddo
             call MPI_SENDRECV(bff_snd(1),nxe-nxs+1,mnpr,ndown,101, &
                               bff_rcv(1),nxe-nxs+1,mnpr,nup  ,101, &
                               ncomw,nstat,nerr)
             do i=nxs,nxe
                ii = i-nxs+1
                p(i,nye+1) = bff_rcv(ii)
             enddo

             do i=nxs,nxe
                ii = i-nxs+1
                bff_snd(ii) = p(i,nye)
             enddo
             call MPI_SENDRECV(bff_snd(1),nxe-nxs+1,mnpr,nup  ,100, &
                               bff_rcv(1),nxe-nxs+1,mnpr,ndown,100, &
                               ncomw,nstat,nerr)
             do i=nxs,nxe
                ii = i-nxs+1
                p(i,nys-1) = bff_rcv(ii)
             enddo

             if(bc == 0)then
                p(nxs-1,nys-1:nye+1) = p(nxe,nys-1:nye+1)
                p(nxe+1,nys-1:nye+1) = p(nxs,nys-1:nye+1)
             else if(bc == -1)then
                p(nxs-1,nys-1:nye+1) = p(nxs+1,nys-1:nye+1)
                p(nxe+1,nys-1:nye+1) = p(nxe-1,nys-1:nye+1)
             else
                write(*,*)'choose bc=0 (periodic) or bc=-1 (reflective)'
                stop
             endif
             !----- end of -----
       
             do j=nys,nye
             do i=nxs,nxe
                ap(i,j) = -p(i,j-1)                    &
                          -p(i-1,j)+f1*p(i,j)-p(i+1,j) &
                          -p(i,j+1)
             enddo
             enddo

             sumr = 0.0
             sum2 = 0.0
             do j=nys,nye
             do i=nxs,nxe
                sumr = sumr+r(i,j)*r(i,j)
                sum2 = sum2+p(i,j)*ap(i,j)
             enddo
             enddo

             bff_snd(1) = sumr
             bff_snd(2) = sum2
             call MPI_ALLREDUCE(bff_snd,bff_rcv,2,mnpr,opsum,ncomw,nerr)
             sumr_g = bff_rcv(1)
             sum2_g = bff_rcv(2)

             av = sumr_g/sum2_g
             
             x(nxs:nxe,nys:nye) = x(nxs:nxe,nys:nye)+av* p(nxs:nxe,nys:nye)
             r(nxs:nxe,nys:nye) = r(nxs:nxe,nys:nye)-av*ap(nxs:nxe,nys:nye)
             
             sum_g = dsqrt(sumr_g)
             if(ite >= ite_max) then
                write(6,*)'********** stop at cgm after ite_max **********'
                stop
             endif
             
             sum1 = 0.0
             do j=nys,nye
             do i=nxs,nxe
                sum1 = sum1+r(i,j)*r(i,j)
             enddo
             enddo
             call MPI_ALLREDUCE(sum1,sum1_g,1,mnpr,opsum,ncomw,nerr)
             bv = sum1_g/sumr_g
             
             p(nxs:nxe,nys:nye) = r(nxs:nxe,nys:nye)+bv*p(nxs:nxe,nys:nye)
             
          enddo
       endif

       gb(l,nxs:nxe,nys:nye) = x(nxs:nxe,nys:nye)

    end do
    
  end subroutine cgm


end module field
