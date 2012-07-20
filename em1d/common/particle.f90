module particle


  implicit none

  private

  public :: particle__init
  public :: particle__solv

  integer, save              :: np, nx, nsp, bc
  real(8), save              :: c, delt
  real(8), save, allocatable :: q(:), r(:)


contains


  subroutine particle__init(npin,nxin,nspin,bcin,qin,rin,cin,deltin)
    integer, intent(in) :: npin, nxin, nspin, bcin
    real(8), intent(in) :: qin(nspin), rin(nspin), cin, deltin

    np   = npin
    nx   = nxin
    nsp  = nspin
    bc   = bcin
    allocate(q(nsp))
    allocate(r(nsp))
    q    = qin
    r    = rin
    c    = cin
    delt = deltin

  end subroutine particle__init


  subroutine particle__solv(gp,up,uf,np2)

    integer, intent(in)  :: np2(1:nx+bc,nsp)
    real(8), intent(in)  :: up(4,np,1:nx+bc,nsp)
    real(8), intent(in)  :: uf(6,0:nx+1)
    real(8), intent(out) :: gp(4,np,1:nx+bc,nsp)
    integer :: i, ii, isp, ih
    real(8) :: pf(6)
    real(8) :: uvm(6)
    real(8) :: dx, dxm
    real(8) :: fac1, fac1r, fac2, fac2r, gam, txxx, bt2

    do isp=1,nsp

       fac1 = q(isp)/r(isp)*0.5*delt
       txxx = fac1*fac1
       fac2 = q(isp)*delt/r(isp)
       do i=1,nx+bc
          do ii=1,np2(i,isp)
             pf(1) = uf(1,i)

             dx = up(1,ii,i,isp)-i
             dxm = 1.-dx
             pf(2) = +dxm*uf(2,i)+dx*uf(2,i+1)
             pf(3) = +dxm*uf(3,i)+dx*uf(3,i+1)
             pf(4) = +dxm*uf(4,i)+dx*uf(4,i+1)

             ih = floor(up(1,ii,i,isp)+0.5)
             dx = up(1,ii,i,isp)+0.5-ih
             dxm = 1.-dx
             pf(5) = +dxm*uf(5,ih-1)+dx*uf(5,ih)
             pf(6) = +dxm*uf(6,ih-1)+dx*uf(6,ih)

             bt2 = pf(1)*pf(1)+pf(2)*pf(2)+pf(3)*pf(3)

             uvm(1) = up(2,ii,i,isp)+fac1*pf(4)
             uvm(2) = up(3,ii,i,isp)+fac1*pf(5)
             uvm(3) = up(4,ii,i,isp)+fac1*pf(6)

             gam = dsqrt(c*c+uvm(1)*uvm(1)+uvm(2)*uvm(2)+uvm(3)*uvm(3))

             fac1r = fac1/gam
             fac2r = fac2/(gam+txxx*bt2/gam)

             uvm(4) = uvm(1)+fac1r*(+uvm(2)*pf(3)-uvm(3)*pf(2))
             uvm(5) = uvm(2)+fac1r*(+uvm(3)*pf(1)-uvm(1)*pf(3))
             uvm(6) = uvm(3)+fac1r*(+uvm(1)*pf(2)-uvm(2)*pf(1))

             uvm(1) = uvm(1)+fac2r*(+uvm(5)*pf(3)-uvm(6)*pf(2))
             uvm(2) = uvm(2)+fac2r*(+uvm(6)*pf(1)-uvm(4)*pf(3))
             uvm(3) = uvm(3)+fac2r*(+uvm(4)*pf(2)-uvm(5)*pf(1))

             gp(2,ii,i,isp) = uvm(1)+fac1*pf(4)
             gp(3,ii,i,isp) = uvm(2)+fac1*pf(5)
             gp(4,ii,i,isp) = uvm(3)+fac1*pf(6)

             gam = dsqrt(1.0+(+gp(2,ii,i,isp)*gp(2,ii,i,isp) &
                              +gp(3,ii,i,isp)*gp(3,ii,i,isp) &
                              +gp(4,ii,i,isp)*gp(4,ii,i,isp))/(c*c))
             gp(1,ii,i,isp) = up(1,ii,i,isp)+gp(2,ii,i,isp)*delt/gam
          enddo
       enddo

    enddo

  end subroutine particle__solv


end module particle
