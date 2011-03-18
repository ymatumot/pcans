program main

  use const
  use boundary
  use fio
  use particle
  use mom_calc

  implicit none

  real(8) :: up(4,np,nxgs:nxge+bc,nsp)
  real(8) :: uf(6,nxgs-1:nxge+1)
  real(8) :: den(nxgs-1:nxge+1,nsp), vel(nxgs-1:nxge+1,3,nsp), temp(nxgs-1:nxge+1,3,nsp)
  integer :: np2(nxgs:nxge+bc,nsp), bcp, it0
  real(8) :: c, q(nsp), r(nsp), delt, delx
  integer :: iwork1, iwork2
  integer :: ndata, idata, isp, irank
  integer :: nxs, nxe
  character(len=64) :: dir
  character(len=64) :: ifile(nproc)

  iwork1 = (nxge-nxgs+1)/nproc
  iwork2 = mod(nxge-nxgs+1,nproc)
  ndata = iargc()
  call getarg(1,dir)
  do idata=2,ndata,nproc

     do irank=0,nproc-1
        nxs = irank*iwork1+nxgs+min(irank,iwork2)
        nxe = nxs+iwork1-1

        if(iwork2 > irank) nxe = nxe+1

        if(irank == nproc-1)then
           if(bc == -1) bcp = -1
           if(bc ==  0) bcp = 0
        else
           bcp = 0
        endif

        call getarg(idata+irank,ifile(irank+1))
        write(*,'(a)')'reading.....  '//trim(dir)//trim(ifile(irank+1))
        call fio__input(up,uf,c,q,r,delt,delx,it0, &
                        np,nsp,np2,nxgs,nxge,nxs,nxe,nproc,bc,bcp, &
                        dir,ifile(irank+1))
     enddo

     call particle__solv(up,uf,                  &
                         c,q,r,0.5*delt,         &
                         np,nsp,np2,nxgs,nxge,bcp)
     call boundary__particle(up,np,nxgs,nxge,nsp,np2,bc)

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
     
     call fio__mom(den,vel,temp,uf,nxgs,nxge,nsp,bc,it0,trim(dir)//'../mom/')
!!$     call fio__psd(up,np,nxgs,nxge,nsp,np2,bc,it0,trim(dir)//'../psd/')
     
  enddo
     
end program main
