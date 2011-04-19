module fio

  implicit none

  private

  public :: fio__input, fio__psd, fio__mom

  logical, save                :: lflag=.true.
  integer, allocatable, public :: np2(:,:)
  integer,              public :: it0, np, nx, nsp, bc
  real(8), allocatable, public :: q(:), r(:), up(:,:,:,:), uf(:,:)
  real(8),              public :: c, delt, delx


contains


  subroutine fio__input(dir,file)

    character(len=*) :: dir, file

    !filename
    open(11,file=trim(dir)//trim(file),form='unformatted')
    read(11)it0,np,nx,nsp,bc,delt,delx,c

    if(lflag)then
       allocate(q(nsp))
       allocate(r(nsp))
       allocate(np2(1:nx+bc,nsp))
       allocate(up(4,np,1:nx+bc,nsp))
       allocate(uf(6,0:nx+1))
       lflag = .false.
    endif

    read(11)np2
    read(11)q
    read(11)r

    !field data
    read(11)uf

    !particle data
    read(11)up

    close(11)

  end subroutine fio__input


  subroutine fio__psd(dir)

    character(len=*), intent(in) :: dir
    integer :: i, ii, isp, ieq
    character(len=256) :: filename

    write(filename,'(a,i6.6,a)')trim(dir),it0,'_'//'psd_i.dat'
    open(90,file=filename,status='unknown')
    write(filename,'(a,i6.6,a)')trim(dir),it0,'_'//'psd_e.dat'
    open(91,file=filename,status='unknown')

    isp=1
    do i=1,nx+bc
       do ii=1,np2(i,isp)
          write(90,'(4e20.5)')(up(ieq,ii,i,isp),ieq=1,4)
       enddo
    enddo

    isp=2
    do i=1,nx+bc
       do ii=1,np2(i,isp)
          write(91,'(4e20.5)')(up(ieq,ii,i,isp),ieq=1,4)
       enddo
    enddo

  end subroutine fio__psd


  subroutine fio__mom(den,vel,temp,dir)

    real(8), intent(inout) :: den(0:nx+1,nsp), vel(0:nx+1,3,nsp), temp(0:nx+1,3,nsp)
    character(len=*), intent(in) :: dir
    integer :: i, j, n_file
    real(8) :: tmp(0:nx+1,6)
    character(len=256) :: filename

    do i=1,nsp
       n_file=10+i-1
       write(filename,'(a,i6.6,a,i1,a)')trim(dir),it0,'_'//'den_',i,'.dat'
       open(n_file,file=filename,status='unknown')
    enddo

    do i=1,nsp
       n_file=10+nsp+(i-1)*3
       write(filename,'(a,i6.6,a,i1,a)')trim(dir),it0,'_'//'Txx_',i,'.dat'
       open(n_file,file=filename,status='unknown')
       n_file=10+nsp+(i-1)*3+1
       write(filename,'(a,i6.6,a,i1,a)')trim(dir),it0,'_'//'Tyy_',i,'.dat'
       open(n_file,file=filename,status='unknown')
       n_file=10+nsp+(i-1)*3+2
       write(filename,'(a,i6.6,a,i1,a)')trim(dir),it0,'_'//'Tzz_',i,'.dat'
       open(n_file,file=filename,status='unknown')
    enddo

    do i=1,nsp
       n_file=10+nsp+3*nsp+(i-1)*3
       write(filename,'(a,i6.6,a,i1,a)')trim(dir),it0,'_'//'vx',i,'.dat'
       open(n_file,file=filename,status='unknown')
       n_file=10+nsp+3*nsp+(i-1)*3+1
       write(filename,'(a,i6.6,a,i1,a)')trim(dir),it0,'_'//'vy',i,'.dat'
       open(n_file,file=filename,status='unknown')
       n_file=10+nsp+3*nsp+(i-1)*3+2
       write(filename,'(a,i6.6,a,i1,a)')trim(dir),it0,'_'//'vz',i,'.dat'
       open(n_file,file=filename,status='unknown')
    enddo

    n_file=10+nsp+3*nsp+3*nsp
    write(filename,'(a,i6.6,a)')trim(dir),it0,'_'//'bx.dat'
    open(n_file,file=filename,status='unknown')
    n_file=10+nsp+3*nsp+3*nsp+1
    write(filename,'(a,i6.6,a)')trim(dir),it0,'_'//'by.dat'
    open(n_file,file=filename,status='unknown')
    n_file=10+nsp+3*nsp+3*nsp+2
    write(filename,'(a,i6.6,a)')trim(dir),it0,'_'//'bz.dat'
    open(n_file,file=filename,status='unknown')

    n_file=10+nsp+3*nsp+3*nsp+3
    write(filename,'(a,i6.6,a)')trim(dir),it0,'_'//'ex.dat'
    open(n_file,file=filename,status='unknown')
    n_file=10+nsp+3*nsp+3*nsp+3+1
    write(filename,'(a,i6.6,a)')trim(dir),it0,'_'//'ey.dat'
    open(n_file,file=filename,status='unknown')
    n_file=10+nsp+3*nsp+3*nsp+3+2
    write(filename,'(a,i6.6,a)')trim(dir),it0,'_'//'ez.dat'
    open(n_file,file=filename,status='unknown')

    !fields at x=i+1/2
    do i=1,nx+bc
       tmp(i,1) = uf(1,i)
       tmp(i,2) = 0.5*(+uf(2,i)+uf(2,i+1))
       tmp(i,3) = 0.5*(+uf(3,i)+uf(3,i+1))
       tmp(i,4) = 0.5*(+uf(4,i)+uf(4,i+1))
       tmp(i,5) = uf(5,i)
       tmp(i,6) = uf(6,i)
    enddo

    !propose vel and temperature
    vel(1:nx+bc,1,1:nsp) = vel(1:nx+bc,1,1:nsp)/den(1:nx+bc,1:nsp)
    vel(1:nx+bc,2,1:nsp) = vel(1:nx+bc,2,1:nsp)/den(1:nx+bc,1:nsp)
    vel(1:nx+bc,3,1:nsp) = vel(1:nx+bc,3,1:nsp)/den(1:nx+bc,1:nsp)
    temp(1:nx+bc,1,1:nsp) = dsqrt(+temp(1:nx+bc,1,1:nsp)/den(1:nx+bc,1:nsp) &
                                  -vel(1:nx+bc,1,1:nsp)*vel(1:nx+bc,1,1:nsp))
    temp(1:nx+bc,2,1:nsp) = dsqrt(+temp(1:nx+bc,2,1:nsp)/den(1:nx+bc,1:nsp) &
                                  -vel(1:nx+bc,2,1:nsp)*vel(1:nx+bc,2,1:nsp))
    temp(1:nx+bc,3,1:nsp) = dsqrt(+temp(1:nx+bc,3,1:nsp)/den(1:nx+bc,1:nsp) &
                                  -vel(1:nx+bc,3,1:nsp)*vel(1:nx+bc,3,1:nsp))

    do j=1,nsp
       n_file=10+j-1
       write(n_file,99)(den(i,j),i=1,nx+bc)
    enddo

    do j=1,nsp
       n_file=10+nsp+(j-1)*3
       write(n_file,99)(temp(i,1,j),i=1,nx+bc)
       n_file=10+nsp+(j-1)*3+1
       write(n_file,99)(temp(i,2,j),i=1,nx+bc)
       n_file=10+nsp+(j-1)*3+2
       write(n_file,99)(temp(i,3,j),i=1,nx+bc)
    enddo

    do j=1,nsp
       n_file=10+nsp+nsp*3+(j-1)*3
       write(n_file,99)(vel(i,1,j),i=1,nx+bc)
       n_file=10+nsp+nsp*3+(j-1)*3+1
       write(n_file,99)(vel(i,2,j),i=1,nx+bc)
       n_file=10+nsp+nsp*3+(j-1)*3+2
       write(n_file,99)(vel(i,3,j),i=1,nx+bc)
    enddo

    do j=1,6
       n_file=10+nsp+nsp*3+nsp*3+(j-1)
       write(n_file,99)(tmp(i,j),i=1,nx+bc)
    enddo

99  format(100000E15.5)

    n_file=10+nsp+nsp*3+nsp*3+5
    do j=10,n_file
       close(j)
    enddo

  end subroutine fio__mom


end module fio
