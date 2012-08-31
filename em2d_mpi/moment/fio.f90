module fio

  implicit none

  private

  public :: fio__input, fio__mom

  logical                      :: lflag=.true.
  integer, public              :: np, nsp, nxgs, nxge, nygs, nyge, nxs, nxe, nys, nye, bc, it0
  integer, public, allocatable :: np2(:,:)
  real(8), public, allocatable :: up(:,:,:,:)
  real(8), public, allocatable :: uf(:,:,:)
  real(8), public, allocatable :: q(:), r(:)
  real(8), public, allocatable :: den(:,:,:), vel(:,:,:,:), temp(:,:,:,:)
  real(8), public              :: c, delt, delx


contains


  subroutine fio__input(nproc,file)
                       
    integer                      :: inproc
    integer, intent(in)          :: nproc
    character(len=*), intent(in) :: file
    real(8), allocatable         :: uftmp(:,:,:)

    open(11,file=trim(file),form='unformatted')

    !parameters
    read(11)it0,np,nxgs,nxge,nygs,nyge,nxs,nxe,nys,nye,nsp,inproc,bc,delt,delx,c
    if(inproc /= nproc)then
       write(*,*)'error in no. of procs'
       stop
    endif

    if(lflag)then
       allocate(q(nsp))
       allocate(r(nsp))
       allocate(uf(6,nxgs-1:nxge+1,nygs-1:nyge+1))
       allocate(den(nxgs-1:nxge+1,nygs-1:nyge+1,nsp))
       allocate(vel(nxgs-1:nxge+1,nygs-1:nyge+1,3,nsp))
       allocate(temp(nxgs-1:nxge+1,nygs-1:nyge+1,6,nsp))

       den(nxgs-1:nxge+1,nygs-1:nyge+1,1:nsp) = 0.0D0
       vel(nxgs-1:nxge+1,nygs-1:nyge+1,1:3,1:nsp) = 0.0D0
       temp(nxgs-1:nxge+1,nygs-1:nyge+1,1:6,1:nsp) = 0.0D0

       lflag = .false.
    endif

    allocate(np2(nys:nye,nsp))
    allocate(up(5,np,nys:nye,nsp))
    allocate(uftmp(6,nxs-1:nxe+1,nys-1:nye+1))

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


  subroutine fio__mom(den,vel,temp,uf,nxgs,nxge,nygs,nyge,nsp,bc,it0,dir)

    integer, intent(in)    :: nxgs, nxge, nygs, nyge, nsp, bc, it0
    real(8), intent(in)    :: uf(6,nxgs-1:nxge+1,nygs-1:nyge+1)
    real(8), intent(inout) :: den(nxgs-1:nxge+1,nygs-1:nyge+1,nsp),   &
                              vel(nxgs-1:nxge+1,nygs-1:nyge+1,3,nsp), &
                              temp(nxgs-1:nxge+1,nygs-1:nyge+1,6,nsp)
    character(len=*), intent(in) :: dir
    integer            :: i, j
    real(8)            :: tmp(nxgs:nxge,nygs:nyge,6)
    character(len=256) :: filename

    write(filename,'(a,i6.6,a)')trim(dir),it0,'_'//'den_i.dat'
    open(10,file=filename,status='unknown')
    write(filename,'(a,i6.6,a)')trim(dir),it0,'_'//'den_e.dat'
    open(11,file=filename,status='unknown')
    write(filename,'(a,i6.6,a)')trim(dir),it0,'_'//'Txx_i.dat'
    open(12,file=filename,status='unknown')
    write(filename,'(a,i6.6,a)')trim(dir),it0,'_'//'Tyy_i.dat'
    open(13,file=filename,status='unknown')
    write(filename,'(a,i6.6,a)')trim(dir),it0,'_'//'Tzz_i.dat'
    open(14,file=filename,status='unknown')
    write(filename,'(a,i6.6,a)')trim(dir),it0,'_'//'Txx_e.dat'
    open(15,file=filename,status='unknown')
    write(filename,'(a,i6.6,a)')trim(dir),it0,'_'//'Tyy_e.dat'
    open(16,file=filename,status='unknown')
    write(filename,'(a,i6.6,a)')trim(dir),it0,'_'//'Tzz_e.dat'
    open(17,file=filename,status='unknown')

    write(filename,'(a,i6.6,a)')trim(dir),it0,'_'//'vxi.dat'
    open(18,file=filename,status='unknown')
    write(filename,'(a,i6.6,a)')trim(dir),it0,'_'//'vyi.dat'
    open(19,file=filename,status='unknown')
    write(filename,'(a,i6.6,a)')trim(dir),it0,'_'//'vzi.dat'
    open(20,file=filename,status='unknown')
    write(filename,'(a,i6.6,a)')trim(dir),it0,'_'//'vxe.dat'
    open(21,file=filename,status='unknown')
    write(filename,'(a,i6.6,a)')trim(dir),it0,'_'//'vye.dat'
    open(22,file=filename,status='unknown')
    write(filename,'(a,i6.6,a)')trim(dir),it0,'_'//'vze.dat'
    open(23,file=filename,status='unknown')

    write(filename,'(a,i6.6,a)')trim(dir),it0,'_'//'bx.dat'
    open(24,file=filename,status='unknown')
    write(filename,'(a,i6.6,a)')trim(dir),it0,'_'//'by.dat'
    open(25,file=filename,status='unknown')
    write(filename,'(a,i6.6,a)')trim(dir),it0,'_'//'bz.dat'
    open(26,file=filename,status='unknown')

    write(filename,'(a,i6.6,a)')trim(dir),it0,'_'//'ex.dat'
    open(27,file=filename,status='unknown')
    write(filename,'(a,i6.6,a)')trim(dir),it0,'_'//'ey.dat'
    open(28,file=filename,status='unknown')
    write(filename,'(a,i6.6,a)')trim(dir),it0,'_'//'ez.dat'
    open(29,file=filename,status='unknown')


    !fields at (i+1/2, j+1/2)
    do j=nygs,nyge
    do i=nxgs,nxge+bc
       tmp(i,j,1) = 0.5*(+uf(1,i,j)+uf(1,i,j+1))
       tmp(i,j,2) = 0.5*(+uf(2,i,j)+uf(2,i+1,j))
       tmp(i,j,3) = 0.25*(+uf(3,i,j)  +uf(3,i+1,j) &
                          +uf(3,i,j+1)+uf(3,i+1,j+1))
       tmp(i,j,4) = 0.5*(+uf(4,i,j)+uf(4,i+1,j))
       tmp(i,j,5) = 0.5*(+uf(5,i,j)+uf(5,i,j+1))
       tmp(i,j,6) = uf(6,i,j)
    enddo
    enddo

    !propose vel and temperature
    vel(nxgs:nxge+bc,nygs:nyge,1,1:nsp) = vel(nxgs:nxge+bc,nygs:nyge,1,1:nsp) &
                                         /den(nxgs:nxge+bc,nygs:nyge,1:nsp)
    vel(nxgs:nxge+bc,nygs:nyge,2,1:nsp) = vel(nxgs:nxge+bc,nygs:nyge,2,1:nsp) &
                                         /den(nxgs:nxge+bc,nygs:nyge,1:nsp)
    vel(nxgs:nxge+bc,nygs:nyge,3,1:nsp) = vel(nxgs:nxge+bc,nygs:nyge,3,1:nsp) &
                                         /den(nxgs:nxge+bc,nygs:nyge,1:nsp)

    temp(nxgs:nxge+bc,nygs:nyge,1,1:nsp) = +temp(nxgs:nxge+bc,nygs:nyge,1,1:nsp) &
                                            /den(nxgs:nxge+bc,nygs:nyge,1:nsp)   &
                                           -vel(nxgs:nxge+bc,nygs:nyge,1,1:nsp)**2
    temp(nxgs:nxge+bc,nygs:nyge,2,1:nsp) = +temp(nxgs:nxge+bc,nygs:nyge,2,1:nsp) &
                                            /den(nxgs:nxge+bc,nygs:nyge,1:nsp)   &
                                           -vel(nxgs:nxge+bc,nygs:nyge,2,1:nsp)**2
    temp(nxgs:nxge+bc,nygs:nyge,3,1:nsp) = +temp(nxgs:nxge+bc,nygs:nyge,3,1:nsp) &
                                            /den(nxgs:nxge+bc,nygs:nyge,1:nsp)   &
                                           -vel(nxgs:nxge+bc,nygs:nyge,3,1:nsp)**2


    do j=nygs,nyge
       write(10,99)(den(i,j,1),i=nxgs,nxge+bc)
       write(11,99)(den(i,j,2),i=nxgs,nxge+bc)
       write(12,99)(temp(i,j,1,1),i=nxgs,nxge+bc)
       write(13,99)(temp(i,j,2,1),i=nxgs,nxge+bc)
       write(14,99)(temp(i,j,3,1),i=nxgs,nxge+bc)
       write(15,99)(temp(i,j,1,2),i=nxgs,nxge+bc)
       write(16,99)(temp(i,j,2,2),i=nxgs,nxge+bc)
       write(17,99)(temp(i,j,3,2),i=nxgs,nxge+bc)
       write(18,99)(vel(i,j,1,1),i=nxgs,nxge+bc)
       write(19,99)(vel(i,j,2,1),i=nxgs,nxge+bc)
       write(20,99)(vel(i,j,3,1),i=nxgs,nxge+bc)
       write(21,99)(vel(i,j,1,2),i=nxgs,nxge+bc)
       write(22,99)(vel(i,j,2,2),i=nxgs,nxge+bc)
       write(23,99)(vel(i,j,3,2),i=nxgs,nxge+bc)
       write(24,99)(tmp(i,j,1),i=nxgs,nxge+bc)
       write(25,99)(tmp(i,j,2),i=nxgs,nxge+bc)
       write(26,99)(tmp(i,j,3),i=nxgs,nxge+bc)
       write(27,99)(tmp(i,j,4),i=nxgs,nxge+bc)
       write(28,99)(tmp(i,j,5),i=nxgs,nxge+bc)
       write(29,99)(tmp(i,j,6),i=nxgs,nxge+bc)
    enddo
    close(10)
    close(11)
    close(12)
    close(13)
    close(14)
    close(15)
    close(16)
    close(17)
    close(18)
    close(19)
    close(20)
    close(21)
    close(22)
    close(23)
    close(24)
    close(25)
    close(26)
    close(27)
    close(28)
    close(29)

99  format(10000E15.5)

  end subroutine fio__mom


end module fio
