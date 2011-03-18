program main

  use const
  use boundary
  use fio
  use particle
  use mom_calc

  implicit none

  integer :: iwork1, iwork2
  integer :: ndata, idata, isp, irank
  integer :: nxs, nxe, nys, nye
  integer :: it0
  integer, allocatable :: np2(:,:) 
  real(8), allocatable :: up(:,:,:,:)
  real(8) :: uf(6,nxgs-1:nxge+1,nygs-1:nyge+1)
  real(8) :: den(nxgs-1:nxge+1,nygs-1:nyge+1,nsp), vel(nxgs-1:nxge+1,nygs-1:nyge+1,3,nsp), &
             temp(nxgs-1:nxge+1,nygs-1:nyge+1,3,nsp)
  real(8) :: c, q(nsp), r(nsp), delt, delx
  character(len=64) :: dir
  character(len=64) :: ifile(nproc)

  iwork1 = (nyge-nygs+1)/nproc
  iwork2 = mod(nyge-nygs+1,nproc)
  ndata = iargc()
  call getarg(1,dir)
  do idata=2,ndata,nproc

     !Initialization
     den(nxgs-1:nxge+1,nygs-1:nyge+1,1:nsp) = 0.0D0
     vel(nxgs-1:nxge+1,nygs-1:nyge+1,1:3,1:nsp) = 0.0D0
     temp(nxgs-1:nxge+1,nygs-1:nyge+1,1:3,1:nsp) = 0.0D0

     do irank=0,nproc-1
        nxs = nxgs
        nxe = nxge
        nys = irank*iwork1+nygs+min(irank,iwork2)
        nye = nys+iwork1-1
        if(iwork2 > irank) nye = nye+1

        allocate(np2(nys:nye,nsp))
        allocate(up(5,np,nys:nye,nsp))

        call getarg(idata+irank,ifile(irank+1))

        write(*,'(a)')'reading.....  '//trim(dir)//trim(ifile(irank+1))

        call fio__input(up,uf,c,q,r,delt,delx,it0, &
                        np,np2,nxgs,nxge,nygs,nyge,nxs,nxe,nys,nye,nsp,nproc,bc, &
                        dir,ifile(irank+1))
        call particle__solv(up,uf,c,q,r,0.5*delt,np,nxgs,nxge,nygs,nyge,nys,nye,nsp,np2)
        call boundary__particle(up,np,nys,nye,nxgs,nxge,nygs,nyge,nsp,np2,bc)
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
     enddo

     call fio__mom(den,vel,temp,uf,nxgs,nxge,nygs,nyge,nsp,bc,it0,trim(dir)//'../mom/')

  enddo

end program main
