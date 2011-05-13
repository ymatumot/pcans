module fio

  implicit none

  private

  public :: fio__input, fio__psd, fio__mom

  logical, save                :: lflag=.true.
  integer, public              :: it0, np, nxgs, nxge, nxs, nxe, nsp, bc, bcp
  integer, public, allocatable :: np2(:,:)
  real(8), public, allocatable :: up(:,:,:,:)
  real(8), public, allocatable :: uf(:,:)
  real(8), public, allocatable :: q(:), r(:)
  real(8), public              :: c, delt, delx


contains


  subroutine fio__input(nproc,dir,file)

    integer, intent(in)          :: nproc
    character(len=*), intent(in) :: dir, file
    integer                      :: inproc, iit0
    integer, allocatable         :: np2tmp(:,:)
    real(8), allocatable         :: uptmp(:,:,:,:),uftmp(:,:)

    open(11,file=trim(dir)//trim(file),form='unformatted')

    !parameters
    read(11)it0,np,nxgs,nxge,nxs,nxe,nsp,inproc,bc,bcp,delt,delx,c

    if(inproc /= nproc)then
       write(*,*)'error in no. of procs'
       stop
    endif

    if(lflag)then
       allocate(q(nsp))
       allocate(r(nsp))
       allocate(np2(nxgs:nxge+bc,nsp))
       allocate(up(4,np,nxgs:nxge+bc,nsp))
       allocate(uf(6,nxgs-1:nxge+1))
       lflag = .false.
    endif

    allocate(np2tmp(nxs:nxe+bcp,nsp))
    allocate(uftmp(6,nxs-1:nxe+1))
    allocate(uptmp(4,np,nxs:nxe+bcp,nsp))

    read(11)np2tmp
    read(11)q
    read(11)r

    !field data
    read(11)uftmp

    !particle data
    read(11)uptmp

    close(11)

    !copy to global variables
    np2(nxs:nxe+bcp,1:nsp) = np2tmp(nxs:nxe+bcp,1:nsp)
    uf(1:6,nxs-1:nxe+1) = uftmp(1:6,nxs-1:nxe+1)
    up(1:4,1:np,nxs:nxe+bcp,1:nsp) = uptmp(1:4,1:np,nxs:nxe+bcp,1:nsp)

    deallocate(np2tmp)
    deallocate(uftmp)
    deallocate(uptmp)

  end subroutine fio__input


  subroutine fio__psd(dir)

    character(len=*), intent(in) :: dir
    integer                      :: i, ii, isp, ieq
    character(len=256)           :: filename

    write(filename,'(a,i6.6,a)')trim(dir),it0,'_'//'psd_i.dat'
    open(90,file=filename,status='unknown')
    write(filename,'(a,i6.6,a)')trim(dir),it0,'_'//'psd_e.dat'
    open(91,file=filename,status='unknown')

    isp=1
    do i=nxgs,nxge+bc
       do ii=1,np2(i,isp)
          write(90,'(4e15.3)')(up(ieq,ii,i,isp),ieq=1,4)
       enddo
    enddo

    isp=2
    do i=nxgs,nxge+bc
       do ii=1,np2(i,isp)
          write(91,'(4e15.3)')(up(ieq,ii,i,isp),ieq=1,4)
       enddo
    enddo

  end subroutine fio__psd


  subroutine fio__mom(den,vel,temp,dir)

    real(8), intent(inout)       :: den(nxgs-1:nxge+1,nsp),   &
                                    vel(nxgs-1:nxge+1,3,nsp), &
                                    temp(nxgs-1:nxge+1,3,nsp)
    character(len=*), intent(in) :: dir
    integer :: i,j,n_file
    real(8) :: tmp(nxgs:nxge,6)
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
    do i=nxgs,nxge+bc
       tmp(i,1) = uf(1,i)
       tmp(i,2) = 0.5*(+uf(2,i)+uf(2,i+1))
       tmp(i,3) = 0.5*(+uf(3,i)+uf(3,i+1))
       tmp(i,4) = 0.5*(+uf(4,i)+uf(4,i+1))
       tmp(i,5) = uf(5,i)
       tmp(i,6) = uf(6,i)
    enddo
    
    !propose vel and temperature
    vel(nxgs:nxge+bc,1,1:nsp) = vel(nxgs:nxge+bc,1,1:nsp)/den(nxgs:nxge+bc,1:nsp)
    vel(nxgs:nxge+bc,2,1:nsp) = vel(nxgs:nxge+bc,2,1:nsp)/den(nxgs:nxge+bc,1:nsp)
    vel(nxgs:nxge+bc,3,1:nsp) = vel(nxgs:nxge+bc,3,1:nsp)/den(nxgs:nxge+bc,1:nsp)
    temp(nxgs:nxge+bc,1,1:nsp) = dsqrt(+temp(nxgs:nxge+bc,1,1:nsp)/den(nxgs:nxge+bc,1:nsp) &
                                        -vel(nxgs:nxge+bc,1,1:nsp)**2)
    temp(nxgs:nxge+bc,2,1:nsp) = dsqrt(+temp(nxgs:nxge+bc,2,1:nsp)/den(nxgs:nxge+bc,1:nsp) &
                                        -vel(nxgs:nxge+bc,2,1:nsp)**2)
    temp(nxgs:nxge+bc,3,1:nsp) = dsqrt(+temp(nxgs:nxge+bc,3,1:nsp)/den(nxgs:nxge+bc,1:nsp) &
                                        -vel(nxgs:nxge+bc,3,1:nsp)**2)

    do j=1,nsp
       n_file=10+j-1
       write(n_file,99)(den(i,j),i=nxgs,nxge+bc)
    enddo

    do j=1,nsp
       n_file=10+nsp+(j-1)*3
       write(n_file,99)(temp(i,1,j),i=nxgs,nxge+bc)
       n_file=10+nsp+(j-1)*3+1
       write(n_file,99)(temp(i,2,j),i=nxgs,nxge+bc)
       n_file=10+nsp+(j-1)*3+2
       write(n_file,99)(temp(i,3,j),i=nxgs,nxge+bc)
    enddo

    do j=1,nsp
       n_file=10+nsp+nsp*3+(j-1)*3
       write(n_file,99)(vel(i,1,j),i=nxgs,nxge+bc)
       n_file=10+nsp+nsp*3+(j-1)*3+1
       write(n_file,99)(vel(i,2,j),i=nxgs,nxge+bc)
       n_file=10+nsp+nsp*3+(j-1)*3+2
       write(n_file,99)(vel(i,3,j),i=nxgs,nxge+bc)
    enddo

    do j=1,6
       n_file=10+nsp+nsp*3+nsp*3+(j-1)
       write(n_file,99)(tmp(i,j),i=nxgs,nxge+bc)
    enddo

99  format(100000E15.5)

    n_file=10+nsp+nsp*3+nsp*3+5
    do j=10,n_file
       close(j)
    enddo

  end subroutine fio__mom


end module fio
