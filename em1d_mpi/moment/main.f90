program main

  use boundary
  use fio
  use particle
  use mom_calc

  implicit none

  logical, save        :: lflag=.true.
  integer              :: nproc, ndata, idata, isp, irank
  integer              :: iargc  !! iargc() - implicit none on Cray
  real(8), allocatable :: den(:,:), vel(:,:,:), temp(:,:,:)
  character(len=64)    :: dir
  character(len=64)    :: ifile

  ndata = iargc()
  call getarg(1,dir)

  write(*,*)'No. of processes?'
  read(*,*)nproc

  do idata=2,ndata,nproc

     do irank=0,nproc-1
        call getarg(idata+irank,ifile)
        write(*,'(a)')'reading.....  '//trim(ifile)

        call fio__input(nproc,ifile)
     enddo

     if(lflag)then
        allocate(den(nxgs-1:nxge+1,nsp))
        allocate(vel(nxgs-1:nxge+1,3,nsp))
        allocate(temp(nxgs-1:nxge+1,3,nsp))
        lflag = .false.
     endif
     
     call particle__solv(up,uf,                  &
                         c,q,r,0.5*delt,         &
                         np,nsp,np2,nxgs,nxge,bcp)

     call mom_calc__den(den,up(1,1:np,nxgs:nxge+bc,1:nsp),np,nxgs,nxge,nsp,np2,bc)
     do isp=1,nsp
        call boundary__den(den(nxgs-1:nxge+1,isp),nxgs,nxge,bc)
     enddo
     
     call mom_calc__vel(vel,up,c,np,nxgs,nxge,nsp,np2,bc)
     do isp=1,nsp
        call boundary__vel(vel(nxgs-1:nxge+1,1:3,isp),nxgs,nxge,bc)
     enddo
     
     call mom_calc__temp(temp,up,c,np,nxgs,nxge,nsp,np2,bc)
     do isp=1,nsp
        call boundary__den(temp(nxgs-1:nxge+1,1,isp),nxgs,nxge,bc)
        call boundary__den(temp(nxgs-1:nxge+1,2,isp),nxgs,nxge,bc)
        call boundary__den(temp(nxgs-1:nxge+1,3,isp),nxgs,nxge,bc)
     enddo
     
     call fio__mom(den,vel,temp,trim(dir)//'/mom/')
     call fio__psd(trim(dir)//'/psd/')
     
  enddo
     
end program main
