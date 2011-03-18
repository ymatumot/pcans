program main

  use const
  use boundary
  use fio
  use particle

  implicit none

  integer :: iwork1, iwork2
  integer :: ndata, idata, irank
  integer :: nxs, nxe, nys, nye
  integer :: it0
  integer, allocatable :: np2(:,:) 
  real(8), allocatable :: up(:,:,:,:)
  real(8) :: uf(6,nxgs-1:nxge+1,nygs-1:nyge+1)
  real(8) :: c, q(nsp), r(nsp), delt, delx, x0, y0, dx, dy
  character(len=64) :: dir, xpos, ypos
  character(len=64) :: ifile(nproc)

  iwork1 = (nyge-nygs+1)/nproc
  iwork2 = mod(nyge-nygs+1,nproc)
  ndata = iargc()
  call getarg(1,xpos)
  call getarg(2,ypos)
  read(xpos,*)x0
  read(ypos,*)y0
  dx = 2.0 !sampling area in the x direction
  dy = 2.0 !sampling area in the y direction
  call getarg(3,dir)

  do idata=4,ndata,nproc
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
        call fio__psd(up,x0,y0,dx*delx,dy*delx,np,nys,nye,nsp,np2,it0,trim(dir)//'../psd/')

        deallocate(np2)
        deallocate(up)
     enddo
  enddo

end program main
