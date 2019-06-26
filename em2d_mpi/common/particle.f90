module particle

  implicit none

  private

  public :: particle__init
  public :: particle__solv
  public :: particle__solv_vay

  logical, save :: is_init = .false.
  integer, save :: np, nsp, nxs, nxe, nys, nye, nsfo
  real(8), save :: c, d_delx, delt
  real(8), allocatable, save :: q(:), r(:)


contains


  subroutine particle__init(np_in,nsp_in,&
                            nxs_in,nxe_in,nys_in,nye_in,nsfo_in, &
                            q_in,r_in,c_in,delx_in,delt_in)

    integer, intent(in) :: np_in, nsp_in
    integer, intent(in) :: nxs_in, nxe_in, nys_in, nye_in, nsfo_in
    real(8), intent(in) :: q_in(nsp_in), r_in(nsp_in), c_in, delx_in, delt_in
 
    np    = np_in
    nsp   = nsp_in
    nxs   = nxs_in
    nxe   = nxe_in
    nys   = nys_in
    nye   = nye_in
    nsfo  = nsfo_in
    c     = c_in
    d_delx = 1./delx_in
    delt  = delt_in

    allocate(q(nsp))
    allocate(r(nsp))
    q     = q_in
    r     = r_in

    is_init = .true.

  end subroutine particle__init


  subroutine particle__solv(gp,up,uf,np2)
                            

    use shape_function, only : sf

    integer, intent(in)  :: np2(nys:nye,nsp)
    real(8), intent(in)  :: up(5,np,nys:nye,nsp)
    real(8), intent(in)  :: uf(6,nxs-2:nxe+2,nys-2:nye+2)
    real(8), intent(out) :: gp(5,np,nys:nye,nsp)
    integer :: i, j, ii, isp, ih, ip, jp
    real(8) :: fac1, fac1r, fac2, fac2r, gam, txxx, bt2
    real(8) :: sh(-2:2,2)
    real(8) :: pf(6)
    real(8) :: uvm(6)
    real(8) :: tmp(1:6,nxs-1:nxe+1,nys-1:nye+1)

    if(.not.is_init)then
       write(6,*)'Initialize first by calling particle__init()'
       stop
    endif

    !fields at (i+1/2, j+1/2)
    do j=nys-1,nye+1
    do i=nxs-1,nxe+1
       tmp(1,i,j) = 0.5*(+uf(1,i,j)+uf(1,i,j+1))
       tmp(2,i,j) = 0.5*(+uf(2,i,j)+uf(2,i+1,j))
       tmp(3,i,j) = 0.25*(+uf(3,i,j)  +uf(3,i+1,j) &
                          +uf(3,i,j+1)+uf(3,i+1,j+1))
       tmp(4,i,j) = 0.5*(+uf(4,i,j)+uf(4,i+1,j))
       tmp(5,i,j) = 0.5*(+uf(5,i,j)+uf(5,i,j+1))
       tmp(6,i,j) = uf(6,i,j)
    enddo
    enddo

    do isp=1,nsp

       fac1 = q(isp)/r(isp)*0.5*delt
       txxx = fac1*fac1
       fac2 = q(isp)/r(isp)*delt
       do j=nys,nye
          do ii=1,np2(j,isp)

             pf(1:6) = 0.D0

             ih = int(up(1,ii,j,isp)*d_delx)

             sh(-2:2,1) = sf(ih,up(1,ii,j,isp)*d_delx-0.5,nsfo)
             sh(-2:2,2) = sf(j ,up(2,ii,j,isp)*d_delx-0.5,nsfo)

             !Bx at (i+1/2, j)
             do jp=-1,1
             do ip=-1,1
                pf(1) = pf(1)+tmp(1,ih+ip,j+jp)*sh(ip,1)*sh(jp,2)
             enddo
             enddo

             !By at (i, j+1/2)
             do jp=-1,1
             do ip=-1,1
                pf(2) = pf(2)+tmp(2,ih+ip,j+jp)*sh(ip,1)*sh(jp,2)
             enddo
             enddo

             !Bz at (i, j)
             do jp=-1,1
             do ip=-1,1
                pf(3) = pf(3)+tmp(3,ih+ip,j+jp)*sh(ip,1)*sh(jp,2)
             enddo
             enddo

             !Ex at (i, j+1/2)
             do jp=-1,1
             do ip=-1,1
                pf(4) = pf(4)+tmp(4,ih+ip,j+jp)*sh(ip,1)*sh(jp,2)
             enddo
             enddo

             !Ey at (i+1/2, j)
             do jp=-1,1
             do ip=-1,1
                pf(5) = pf(5)+tmp(5,ih+ip,j+jp)*sh(ip,1)*sh(jp,2)
             enddo
             enddo

             !Ez at (i+1/2, j+1/2)
             do jp=-1,1
             do ip=-1,1
                pf(6) = pf(6)+tmp(6,ih+ip,j+jp)*sh(ip,1)*sh(jp,2)
             enddo
             enddo

             bt2 = pf(1)*pf(1)+pf(2)*pf(2)+pf(3)*pf(3)

             uvm(1) = up(3,ii,j,isp)+fac1*pf(4)
             uvm(2) = up(4,ii,j,isp)+fac1*pf(5)
             uvm(3) = up(5,ii,j,isp)+fac1*pf(6)

             gam = sqrt(c*c+uvm(1)*uvm(1)+uvm(2)*uvm(2)+uvm(3)*uvm(3))
             fac1r = fac1/gam
             fac2r = fac2/(gam+txxx*bt2/gam)

             uvm(4) = uvm(1)+fac1r*(+uvm(2)*pf(3)-uvm(3)*pf(2))
             uvm(5) = uvm(2)+fac1r*(+uvm(3)*pf(1)-uvm(1)*pf(3))
             uvm(6) = uvm(3)+fac1r*(+uvm(1)*pf(2)-uvm(2)*pf(1))
    
             uvm(1) = uvm(1)+fac2r*(+uvm(5)*pf(3)-uvm(6)*pf(2))
             uvm(2) = uvm(2)+fac2r*(+uvm(6)*pf(1)-uvm(4)*pf(3))
             uvm(3) = uvm(3)+fac2r*(+uvm(4)*pf(2)-uvm(5)*pf(1))

             gp(3,ii,j,isp) = uvm(1)+fac1*pf(4)
             gp(4,ii,j,isp) = uvm(2)+fac1*pf(5)
             gp(5,ii,j,isp) = uvm(3)+fac1*pf(6)

             gam = 1./sqrt(1.0D0+(+gp(3,ii,j,isp)*gp(3,ii,j,isp) &
                                  +gp(4,ii,j,isp)*gp(4,ii,j,isp) &
                                  +gp(5,ii,j,isp)*gp(5,ii,j,isp))/(c*c))
             gp(1,ii,j,isp) = up(1,ii,j,isp)+gp(3,ii,j,isp)*delt*gam
             gp(2,ii,j,isp) = up(2,ii,j,isp)+gp(4,ii,j,isp)*delt*gam
          enddo
       enddo
    enddo

  end subroutine particle__solv


! Vay solver [PoP 15, 056701 (2008)] written by SZ
  subroutine particle__solv_vay(gp,up,uf,np2)

    use shape_function, only : sf

    integer, intent(in)  :: np2(nys:nye,nsp)
    real(8), intent(in)  :: up(5,np,nys:nye,nsp)
    real(8), intent(in)  :: uf(6,nxs-2:nxe+2,nys-2:nye+2)
    real(8), intent(out) :: gp(5,np,nys:nye,nsp)
    integer :: i, j, ii, isp, ih, ip, jp
    real(8) :: fac1, fac1r, fac2, txxx, gam, gam2, tau2, ua, sigma
    real(8) :: sh(-2:2,2)
    real(8) :: pf(6)
    real(8) :: uvm(6)
    real(8) :: tmp(1:6,nxs-1:nxe+1,nys-1:nye+1)

    if(.not.is_init)then
       write(6,*)'Initialize first by calling particle__init()'
       stop
    endif

    !fields at (i+1/2, j+1/2)
    do j=nys-1,nye+1
    do i=nxs-1,nxe+1
       tmp(1,i,j) = 0.5*(+uf(1,i,j)+uf(1,i,j+1))
       tmp(2,i,j) = 0.5*(+uf(2,i,j)+uf(2,i+1,j))
       tmp(3,i,j) = 0.25*(+uf(3,i,j)  +uf(3,i+1,j) &
                          +uf(3,i,j+1)+uf(3,i+1,j+1))
       tmp(4,i,j) = 0.5*(+uf(4,i,j)+uf(4,i+1,j))
       tmp(5,i,j) = 0.5*(+uf(5,i,j)+uf(5,i,j+1))
       tmp(6,i,j) = uf(6,i,j)
    enddo
    enddo

    do isp=1,nsp

       fac1 = q(isp)/r(isp)*0.5*delt
       fac2 = q(isp)/r(isp)*delt
       do j=nys,nye
          do ii=1,np2(j,isp)

             pf(1:6) = 0.D0

             ih = int(up(1,ii,j,isp)*d_delx)

             sh(-2:2,1) = sf(ih,up(1,ii,j,isp)*d_delx-0.5,nsfo)
             sh(-2:2,2) = sf(j ,up(2,ii,j,isp)*d_delx-0.5,nsfo)

             !Bx at (i+1/2, j)
             do jp=-1,1
             do ip=-1,1
                pf(1) = pf(1)+tmp(1,ih+ip,j+jp)*sh(ip,1)*sh(jp,2)
             enddo
             enddo

             !By at (i, j+1/2)
             do jp=-1,1
             do ip=-1,1
                pf(2) = pf(2)+tmp(2,ih+ip,j+jp)*sh(ip,1)*sh(jp,2)
             enddo
             enddo

             !Bz at (i, j)
             do jp=-1,1
             do ip=-1,1
                pf(3) = pf(3)+tmp(3,ih+ip,j+jp)*sh(ip,1)*sh(jp,2)
             enddo
             enddo

             !Ex at (i, j+1/2)
             do jp=-1,1
             do ip=-1,1
                pf(4) = pf(4)+tmp(4,ih+ip,j+jp)*sh(ip,1)*sh(jp,2)
             enddo
             enddo

             !Ey at (i+1/2, j)
             do jp=-1,1
             do ip=-1,1
                pf(5) = pf(5)+tmp(5,ih+ip,j+jp)*sh(ip,1)*sh(jp,2)
             enddo
             enddo

             !Ez at (i+1/2, j+1/2)
             do jp=-1,1
             do ip=-1,1
                pf(6) = pf(6)+tmp(6,ih+ip,j+jp)*sh(ip,1)*sh(jp,2)
             enddo
             enddo

             uvm(1) = up(3,ii,j,isp)
             uvm(2) = up(4,ii,j,isp)
             uvm(3) = up(5,ii,j,isp)

             gam = sqrt(c*c+uvm(1)*uvm(1)+uvm(2)*uvm(2)+uvm(3)*uvm(3))
             fac1r = fac1/gam
             uvm(4) = uvm(1) + fac2*pf(4) + fac1r*(+uvm(2)*pf(3)-uvm(3)*pf(2))
             uvm(5) = uvm(2) + fac2*pf(5) + fac1r*(+uvm(3)*pf(1)-uvm(1)*pf(3))
             uvm(6) = uvm(3) + fac2*pf(6) + fac1r*(+uvm(1)*pf(2)-uvm(2)*pf(1))

             pf(1:3) = fac1*pf(1:3)/c  ! tau vector
             tau2  = dot_product(pf(1:3),pf(1:3))
             ua    = dot_product(uvm(4:6),pf(1:3))/c  ! u* = u'.tau/c
             sigma = 1.d0 + dot_product(uvm(4:6),uvm(4:6))/(c*c)-tau2
             gam2  = 0.5d0*(sigma+sqrt(sigma*sigma+4.0*(tau2+ua*ua)))
             gam   = sqrt(gam2)

             txxx  = 1.d0/(tau2+gam2)  ! s
             gp(3,ii,j,isp) = txxx*(gam2*uvm(4)+c*ua*pf(1)+gam*(uvm(5)*pf(3)-uvm(6)*pf(2)))
             gp(4,ii,j,isp) = txxx*(gam2*uvm(5)+c*ua*pf(2)+gam*(uvm(6)*pf(1)-uvm(4)*pf(3)))
             gp(5,ii,j,isp) = txxx*(gam2*uvm(6)+c*ua*pf(3)+gam*(uvm(4)*pf(2)-uvm(5)*pf(1)))

             gam = 1./sqrt(1.0D0+(+gp(3,ii,j,isp)*gp(3,ii,j,isp) &
                                  +gp(4,ii,j,isp)*gp(4,ii,j,isp) &
                                  +gp(5,ii,j,isp)*gp(5,ii,j,isp))/(c*c))
             gp(1,ii,j,isp) = up(1,ii,j,isp)+gp(3,ii,j,isp)*delt*gam
             gp(2,ii,j,isp) = up(2,ii,j,isp)+gp(4,ii,j,isp)*delt*gam
          enddo
       enddo
    enddo

  end subroutine particle__solv_vay


end module particle
