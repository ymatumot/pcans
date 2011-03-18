module fio

  implicit none

  private

  public :: fio__input, fio__psd


contains


  subroutine fio__input(up,uf,c,q,r,delt,delx,it0,                               &
                        np,np2,nxgs,nxge,nygs,nyge,nxs,nxe,nys,nye,nsp,nproc,bc, &
                        dir,file)
                       
    integer, intent(in)  :: np, nsp, nxgs, nxge, nygs, nyge, nxs, nxe, nys, nye, nproc, bc
    character(len=*), intent(in) :: dir, file
    integer, intent(out) :: np2(nys:nye,nsp), it0
    real(8), intent(out) :: up(5,np,nys:nye,nsp)
    real(8), intent(out) :: uf(6,nxgs-1:nxge+1,nygs-1:nyge+1)
    real(8), intent(out) :: c, q(nsp), r(nsp), delt, delx
    integer              :: inp, inxgs, inxge, inygs, inyge, insp, inproc, ibc
    real(8), allocatable :: uftmp(:,:,:)

    allocate(uftmp(6,nxs-1:nxe+1,nys-1:nye+1))

    open(11,file=trim(dir)//trim(file),form='unformatted')

    !parameters
    read(11)it0,inp,inxgs,inxge,inygs,inyge,insp,inproc,ibc,delt,delx,c
    if((inxgs /= nxgs) .or. (inxge /= nxge) .or. (inygs /= nygs) .or. (inyge /= nyge) &
       .or. (inproc /= nproc) .or. (inp /= np) .or. (insp /= nsp) .or. (ibc /= bc))then
       write(6,*) '** parameter mismatch **'
       stop
    endif
    read(11)np2
    read(11)q
    read(11)r

    !field data
    read(11)uftmp
    !particle data
    read(11)up

    close(11)

    !copy to global variables
    uf(1:6,nxs-1:nxe+1,nys-1:nye+1) = uftmp(1:6,nxs-1:nxe+1,nys-1:nye+1)

    deallocate(uftmp)

  end subroutine fio__input


  subroutine fio__psd(up,x0,y0,dx,dy,np,nys,nye,nsp,np2,it0,dir)

    integer, intent(in) :: np, nys, nye, nsp, np2(nys:nye,nsp), it0
    real(8), intent(in) :: up(5,np,nys:nye,nsp)
    real(8), intent(in) :: x0, y0, dx, dy
    character(len=*), intent(in) :: dir
    integer            :: j, ii, isp, ieq
    character(len=256) :: filename

    write(filename,'(a,i6.6,a,i4.4,a,i4.4,a)')trim(dir),it0,'_',floor(x0),'-',floor(y0),'_'//'psd_i.dat'
    open(90,file=filename,status='unknown')
    write(filename,'(a,i6.6,a,i4.4,a,i4.4,a)')trim(dir),it0,'_',floor(x0),'-',floor(y0),'_'//'psd_e.dat'
    open(91,file=filename,status='unknown')

    isp=1
    do j=nys,nye
       do ii=1,np2(j,isp)
          if(up(1,ii,j,isp) >= x0-dx .and. up(1,ii,j,isp) <= x0+dx .and. &
             up(2,ii,j,isp) >= y0-dy .and. up(2,ii,j,isp) <= y0+dy)then
             write(90,'(5e17.7)')(up(ieq,ii,j,isp),ieq=1,5)
          endif
       enddo
    enddo

    isp=2
    do j=nys,nye
       do ii=1,np2(j,isp)
          if(up(1,ii,j,isp) >= x0-dx .and. up(1,ii,j,isp) <= x0+dx .and. &
             up(2,ii,j,isp) >= y0-dy .and. up(2,ii,j,isp) <= y0+dy)then
             write(91,'(5e17.7)')(up(ieq,ii,j,isp),ieq=1,5)
          endif
       enddo
    enddo

  end subroutine fio__psd


end module fio
