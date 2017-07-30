program main

  use boundary
  use fio
  use particle
  use mom_calc

  implicit none

  integer           :: nproc, ndata, idata, isp, irank
  integer           :: iargc  !! iargc() - implicit none on Cray
  character(len=64) :: dir
  character(len=64) :: ifile

  ndata = iargc()
  call getarg(1,dir)

  write(*,*)'No. of processes?'
  read(*,*)nproc

  do idata=2,ndata,nproc

     do irank=0,nproc-1
        call getarg(idata+irank,ifile)
        write(*,'(a)')'reading.....  '//trim(ifile)

        call fio__input(nproc,ifile)

        call particle__solv(up,uf,c,q,r,0.5*delt,np,nxgs,nxge,nygs,nyge,nys,nye,nsp,np2)

        call mom_calc__den(den,up,np,nys,nye,nxgs,nxge,nygs,nyge,nsp,np2)
        call mom_calc__vel(vel,up,c,np,nys,nye,nxgs,nxge,nygs,nyge,nsp,np2)
        call mom_calc__temp(temp,up,c,np,nys,nye,nxgs,nxge,nygs,nyge,nsp,np2)

        deallocate(np2)
        deallocate(up)
     enddo

     do isp=1,nsp
        call boundary__den(den(nxgs-1:nxge+1,nygs-1:nyge+1,isp),nxgs,nxge,nygs,nyge,bc)
     enddo
     do isp=1,nsp
        call boundary__vel(vel(nxgs-1:nxge+1,nygs-1:nyge+1,1:3,isp),nxgs,nxge,nygs,nyge,bc)
     enddo
     do isp=1,nsp
        call boundary__den(temp(nxgs-1:nxge+1,nygs-1:nyge+1,1,isp),nxgs,nxge,nygs,nyge,bc)
        call boundary__den(temp(nxgs-1:nxge+1,nygs-1:nyge+1,2,isp),nxgs,nxge,nygs,nyge,bc)
        call boundary__den(temp(nxgs-1:nxge+1,nygs-1:nyge+1,3,isp),nxgs,nxge,nygs,nyge,bc)
        call boundary__den(temp(nxgs-1:nxge+1,nygs-1:nyge+1,4,isp),nxgs,nxge,nygs,nyge,bc)
        call boundary__den(temp(nxgs-1:nxge+1,nygs-1:nyge+1,5,isp),nxgs,nxge,nygs,nyge,bc)
        call boundary__den(temp(nxgs-1:nxge+1,nygs-1:nyge+1,6,isp),nxgs,nxge,nygs,nyge,bc)
     enddo

     call fio__mom(den,vel,temp,uf,nxgs,nxge,nygs,nyge,nsp,bc,it0,trim(dir)//'/mom/')

     !memory clear
     den(nxgs-1:nxge+1,nygs-1:nyge+1,1:nsp) = 0.0D0
     vel(nxgs-1:nxge+1,nygs-1:nyge+1,1:3,1:nsp) = 0.0D0
     temp(nxgs-1:nxge+1,nygs-1:nyge+1,1:6,1:nsp) = 0.0D0

  enddo

end program main
