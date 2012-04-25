program main

  use boundary
  use fio
  use particle
  use mom_calc

  implicit none

  logical, save        :: lflag=.true.
  integer              :: ndata, idata, isp
  real(8), allocatable :: den(:,:), vel(:,:,:), temp(:,:,:)
  character(len=64)    :: dir
  character(len=64)    :: ifile

  ndata = iargc()
  call getarg(1,dir)
  
  do idata=2,ndata

     call getarg(idata,ifile)
     write(*,'(a)')'reading....  '//trim(dir)//trim(ifile)

     call fio__input(dir,ifile)

     if(lflag)then
        allocate(den(0:nx+1,nsp))
        allocate(vel(0:nx+1,3,nsp))
        allocate(temp(0:nx+1,3,nsp))
        lflag = .false.
     endif
     
     call particle__solv(up,uf,c,q,r,0.5*delt,np2,np,nx,nsp,bc)

     call mom_calc__den(den,up(1,1:np,1:nx+bc,1:nsp),np,nx,nsp,np2,bc)
     do isp=1,nsp
        call boundary__den(den(0:nx+1,isp),nx,bc)
     enddo

     call mom_calc__vel(vel,up,c,np,nx,nsp,np2,bc)
     do isp=1,nsp
        call boundary__vel(vel(0:nx+1,1:3,isp),nx,bc)
     enddo

     call mom_calc__temp(temp,up,c,np,nx,nsp,np2,bc)
     do isp=1,nsp
        call boundary__den(temp(0:nx+1,1,isp),nx,bc)
        call boundary__den(temp(0:nx+1,2,isp),nx,bc)
        call boundary__den(temp(0:nx+1,3,isp),nx,bc)
     enddo

     call fio__mom(den,vel,temp,trim(dir)//'../mom/')
     call fio__psd(trim(dir)//'../psd/')
  enddo

end program main
