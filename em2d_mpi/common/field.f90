module field

  implicit none

  private

  public :: field__init
  public :: field__fdtd_i

  logical, save :: is_init = .false.
  integer, save :: np, nsp, nxs, nxe, nys, nye, nsfo, bc
  integer, save :: nup, ndown, mnpr, opsum, ncomw
  integer       :: nerr
  integer, allocatable :: nstat(:)
  real(8), parameter :: pi = 4.0*atan(1.0d0)
  real(8), save :: c, delx, delt, gfac, d_delx, d_delt
  real(8), allocatable, save :: q(:)


contains

  
  subroutine field__init(np_in,nsp_in,&
                         nxs_in,nxe_in,nys_in,nye_in,nsfo_in,bc_in, &
                         q_in,c_in,delx_in,delt_in,gfac_in,         &
                         nup_in,ndown_in,mnpr_in,opsum_in,ncomw_in,nerr_in,nstat_in)

    integer, intent(in) :: np_in, nsp_in
    integer, intent(in) :: nxs_in, nxe_in, nys_in, nye_in, nsfo_in, bc_in
    integer, intent(in) :: nup_in, ndown_in, mnpr_in, opsum_in, ncomw_in
    integer, intent(in) :: nerr_in, nstat_in(:)
    real(8), intent(in) :: q_in(nsp_in), c_in, delx_in, delt_in, gfac_in
 
    np    = np_in
    nsp   = nsp_in
    nxs   = nxs_in
    nxe   = nxe_in
    nys   = nys_in
    nye   = nye_in
    nsfo  = nsfo_in
    bc    = bc_in
    nup   = nup_in
    ndown = ndown_in
    mnpr  = mnpr_in
    opsum = opsum_in
    ncomw = ncomw_in
    nerr  = nerr_in
    c     = c_in
    delx  = delx_in
    delt  = delt_in
    gfac  = gfac_in
    d_delx = 1./delx
    d_delt = 1./delt

    allocate(nstat(size(nstat_in)))
    nstat = nstat_in
    allocate(q(nsp))
    q     = q_in

    is_init = .true.

  end subroutine field__init


  subroutine field__fdtd_i(uf,up,gp,np2)

    use boundary, only : boundary__dfield, boundary__curre
 
    integer, intent(in)    :: np2(nys:nye,nsp)
    real(8), intent(in)    :: gp(5,np,nys:nye,nsp)
    real(8), intent(inout) :: up(5,np,nys:nye,nsp)
    real(8), intent(inout) :: uf(6,nxs-2:nxe+2,nys-2:nye+2)
    logical, save              :: lflag=.true.
    integer                    :: ii, i, j, isp
    real(8)                    :: f1, f2, f3
    real(8)                    :: uj(3,nxs-2:nxe+2,nys-2:nye+2)
    real(8), save, allocatable :: df(:,:,:), gkl(:,:,:)

    if(.not.is_init)then
       write(6,*)'Initialize first by calling field__init()'
       stop
    endif

    if(lflag)then
       allocate(df(6,nxs-2:nxe+2,nys-2:nye+2))
       allocate(gkl(3,nxs:nxe,nys:nye))
       df(1:6,nxs-2:nxe+2,nys-2:nye+2) = 0.0D0
       gkl(1:3,nxs:nxe,nys:nye) = 0.0D0
       lflag = .false.
    endif

    call ele_cur(uj,up,gp,np2)
    call boundary__curre(uj)

    !calculation
    f1 = c*delt*d_delx
    f2 = gfac*f1*f1
    f3 = 4.0*pi*delx/c
    do j=nys,nye
    do i=nxs,nxe
       gkl(1,i,j) = +f2*(+uf(1,i,j-1)                          &
                         +uf(1,i-1,j)-4.*uf(1,i,j)+uf(1,i+1,j) &
                         +uf(1,i,j+1)                          &
                         +f3*(-uj(3,i,j-1)+uj(3,i,j)) )        &
                    -f1*(-uf(6,i,j-1)+uf(6,i,j))
       gkl(2,i,j) = +f2*(+uf(2,i,j-1)                          &
                         +uf(2,i-1,j)-4.*uf(2,i,j)+uf(2,i+1,j) &
                         +uf(2,i,j+1)                          &
                         -f3*(-uj(3,i-1,j)+uj(3,i,j)) )        &
                    +f1*(-uf(6,i-1,j)+uf(6,i,j))
       gkl(3,i,j) = +f2*(+uf(3,i,j-1)                          &
                         +uf(3,i-1,j)-4.*uf(3,i,j)+uf(3,i+1,j) &
                         +uf(3,i,j+1)                          &
                         +f3*(-uj(2,i-1,j)+uj(2,i,j)           &
                              +uj(1,i,j-1)-uj(1,i,j)) )        &
                    -f1*(-uf(5,i-1,j)+uf(5,i,j)+uf(4,i,j-1)-uf(4,i,j))
    enddo
    enddo

    !solve  < bx, by & bz >
    call cgm(df,gkl)
    call boundary__dfield(df)

    !solve  < ex, ey & ez >
    do j=nys,nye
    do i=nxs,nxe
       df(4,i,j) = +f1*(+gfac*(-df(3,i,j)+df(3,i,j+1))   &
                        +     (-uf(3,i,j)+uf(3,i,j+1)) ) &
                   -4.*pi*delt*uj(1,i,j)
       df(5,i,j) = -f1*(+gfac*(-df(3,i,j)+df(3,i+1,j))   &
                        +     (-uf(3,i,j)+uf(3,i+1,j)) ) &
                   -4.*pi*delt*uj(2,i,j)

       df(6,i,j) = +f1*(+gfac*(-df(2,i,j)+df(2,i+1,j)    &
                               +df(1,i,j)-df(1,i,j+1))   &
                        +     (-uf(2,i,j)+uf(2,i+1,j)    &
                               +uf(1,i,j)-uf(1,i,j+1)) ) &
                   -4.*pi*delt*uj(3,i,j)
    enddo
    enddo

    call boundary__dfield(df)

    !===== Update fields and particles ======
    uf = uf + df

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


  subroutine ele_cur(uj,up,gp,np2)

    use shape_function, only : sf

    integer, intent(in)  :: np2(nys:nye,nsp)
    real(8), intent(in)  :: up(5,np,nys:nye,nsp), gp(5,np,nys:nye,nsp)
    real(8), intent(out) :: uj(3,nxs-2:nxe+2,nys-2:nye+2)
    
    integer            :: ii, j, isp, i1, j1, i2, j2, iinc, jinc, ip, jp
    real(8), parameter :: fac = 1.D0/3.D0
    real(8)            :: gamp
    real(8)            :: s0(-2:2,2), s1(-2:2,2), ds(-2:2,2), ujp(3,-2:3,-2:3)

    uj(1:3,nxs-2:nxe+2,nys-2:nye+2) = 0.D0

    !--------------Charge Conservation Method -------------!
    !---- Density Decomposition (Esirkepov, CPC, 2001) ----!
    do isp=1,nsp
       do j=nys,nye
          do ii=1,np2(j,isp)

             ujp(1:3,-2:3,-2:3) = 0.D0

             i1 = int(up(1,ii,j,isp)*d_delx)
             i2 = int(gp(1,ii,j,isp)*d_delx)
             j2 = int(gp(2,ii,j,isp)*d_delx)
             iinc = i2-i1
             jinc = j2-j

             gamp = 1./sqrt(1.D0+(+gp(3,ii,j,isp)*gp(3,ii,j,isp) &
                                  +gp(4,ii,j,isp)*gp(4,ii,j,isp) &
                                  +gp(5,ii,j,isp)*gp(5,ii,j,isp))/(c*c) )

             s0(-2:2,1) = sf(i1,up(1,ii,j,isp)*d_delx-0.5,nsfo)
             s0(-2:2,2) = sf(j ,up(2,ii,j,isp)*d_delx-0.5,nsfo)

             s1(-2:2,1) = sf(i2,gp(1,ii,j,isp)*d_delx-0.5,nsfo)
             s1(-2:2,2) = sf(j2,gp(2,ii,j,isp)*d_delx-0.5,nsfo)

             ds(-2:2,1) = cshift(s1(-2:2,1),-iinc,1)-s0(-2:2,1)
             ds(-2:2,2) = cshift(s1(-2:2,2),-jinc,1)-s0(-2:2,2)

             do jp=-1+min(jinc,0),1+max(jinc,0)
             do ip=-1+min(iinc,0),1+max(iinc,0)
                ujp(1,ip+1,jp) = ujp(1,ip,jp) &
                                -q(isp)*delx*d_delt*ds(ip,1)*(s0(jp,2)+0.5*ds(jp,2))

                ujp(2,ip,jp+1) = ujp(2,ip,jp) &
                                -q(isp)*delx*d_delt*ds(jp,2)*(s0(ip,1)+0.5*ds(ip,1))

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


  subroutine cgm(df,gkl)

    use boundary, only : boundary__phi

    !-----------------------------------------------------------------------
    !  #  conjugate gradient method 
    !  #  this routine will be stopped after iteration number = ite_max
    !-----------------------------------------------------------------------
    real(8), intent(in)    :: gkl(3,nxs:nxe,nys:nye)
    real(8), intent(inout) :: df(6,nxs-2:nxe+2,nys-2:nye+2)
    integer, parameter :: ite_max = 100 ! maximum number of iteration
    integer            :: i, ii, j, l, ite, bcc
    real(8), parameter :: err = 1D-6 
    real(8)            :: f1, f2, eps, sumr, sum, sum1, sum2, av, bv
    real(8)            :: sumr_g, sum_g, sum1_g, sum2_g
    real(8)            :: phi(nxs-1:nxe+1,nys-1:nye+1), p(nxs-1:nxe+1,nys-1:nye+1)
    real(8)            :: r(nxs:nxe,nys:nye), b(nxs:nxe,nys:nye)
    real(8)            :: ap(nxs:nxe,nys:nye)
    real(8)            :: bff_snd(2), bff_rcv(2)

    f1 = 4.0D0+(delx/(c*delt*gfac))**2
    f2 = (delx/(c*delt*gfac))**2

    do l=1,3

       select case(l)
         case(1)
           bcc = bc
         case(2,3)
           bcc = 0
       endselect

       ! initial guess
       ite = 0
       sum = 0.0D0
       do j=nys,nye
       do i=nxs,nxe+bcc
          phi(i,j) = df(l,i,j)
          b(i,j) = f2*gkl(l,i,j)
          sum = sum+b(i,j)*b(i,j)
       enddo
       enddo

       call MPI_ALLREDUCE(sum,sum_g,1,mnpr,opsum,ncomw,nerr)

       eps = sqrt(sum_g)*err

       !------ boundary condition of phi ------
       call boundary__phi(phi,l)
       !------ end of -----

       sumr = 0.0D0
       do j=nys,nye
       do i=nxs,nxe+bcc
          r(i,j) = b(i,j)+phi(i,j-1)                    &
                         +phi(i-1,j)-f1*phi(i,j)+phi(i+1,j) &
                         +phi(i,j+1)
          p(i,j) = r(i,j)
          sumr = sumr+r(i,j)*r(i,j)
       enddo
       enddo

       call MPI_ALLREDUCE(sumr,sumr_g,1,mnpr,opsum,ncomw,nerr)

       if(sqrt(sumr_g) > eps)then
       
          do while(sum_g > eps)
             
             ite = ite+1

             !------boundary condition of p------
             call boundary__phi(p,l)
             !----- end of -----
       
             sumr = 0.0D0
             sum2 = 0.0D0
             do j=nys,nye
             do i=nxs,nxe+bcc
                ap(i,j) = -p(i,j-1)                    &
                          -p(i-1,j)+f1*p(i,j)-p(i+1,j) &
                          -p(i,j+1)
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
             
             phi(nxs:nxe+bcc,nys:nye) = phi(nxs:nxe+bcc,nys:nye)+av* p(nxs:nxe+bcc,nys:nye)
             r(nxs:nxe+bcc,nys:nye) = r(nxs:nxe+bcc,nys:nye)-av*ap(nxs:nxe+bcc,nys:nye)
             
             sum_g = sqrt(sumr_g)
             if(ite >= ite_max) then
                write(6,*)'********** stop at cgm after ite_max **********'
                stop
             endif
             
             sum1 = 0.0
             do j=nys,nye
             do i=nxs,nxe+bcc
                sum1 = sum1+r(i,j)*r(i,j)
             enddo
             enddo
             call MPI_ALLREDUCE(sum1,sum1_g,1,mnpr,opsum,ncomw,nerr)
             bv = sum1_g/sumr_g
             
             p(nxs:nxe+bcc,nys:nye) = r(nxs:nxe+bcc,nys:nye)+bv*p(nxs:nxe+bcc,nys:nye)
             
          enddo
       endif

       df(l,nxs:nxe+bcc,nys:nye) = phi(nxs:nxe+bcc,nys:nye)

    end do

    
  end subroutine cgm


end module field
