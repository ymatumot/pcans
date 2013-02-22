program main

  use const
  use mpi_set
  use init
  use boundary
  use fio
  use particle
  use field

  implicit none

  integer :: it=0

!**********************************************************************c
!
!    two-dimensional electromagnetic plasma simulation code
!
!    written by M Hoshino,  ISAS, 1984/09/12
!    revised  1985/03/08  1985/04/05  1997/05/06
!    revised for CANS    (by Y. Matsumoto, STEL)  2004/06/22
!    re-written in F90   (by Y. Matsumoto, STEL)  2008/10/21
!    MPI parallelization (by Y. Matsumoto, STEL)  2009/4/1
!    2-D code            (by Y. Matsumoto, STEL)  2009/6/5
!
!**********************************************************************c

  call init__set_param
  call fio__energy(up,uf,                         &
                   np,nsp,np2,nxs,nxe,nys,nye,bc, &
                   c,r,delt,0,it0,dir,file12,     &
                   nroot,nrank,mnpr,opsum,ncomw,nerr)
  call fio__output(up,uf,np,nxgs,nxge,nygs,nyge,nxs,nxe,nys,nye,nsp,np2,bc,nproc,nrank, &
                   c,q,r,delt,delx,0,it0,dir)

  do it=1,itmax-it0

     call particle__solv(gp,up,uf,                        &
                         np,nsp,np2,nxs,nxe,nys,nye,nsfo, &
                         c,q,r,delt,delx)                 
     call field__fdtd_i(uf,up,gp,                           &
                        np,nsp,np2,nxs,nxe,nys,nye,nsfo,bc, &
                        q,c,delx,delt,gfac,                 &
                        nup,ndown,mnpr,opsum,nstat,ncomw,nerr)
     call boundary__particle(up,                                        &
                             np,nsp,np2,nxgs,nxge,nygs,nyge,nys,nye,bc, &
                             nup,ndown,nstat,mnpi,mnpr,ncomw,nerr)

     if (mod(it+it0,intvl1) == 0) then
        call fio__output(up,uf,np,nxgs,nxge,nygs,nyge,nxs,nxe,nys,nye,nsp,np2,bc,nproc,nrank, &
             c,q,r,delt,delx,it,it0,dir)
     endif
     if (mod(it+it0,intvl2) == 0) then
        call fio__energy(up,uf,                         &
                         np,nsp,np2,nxs,nxe,nys,nye,bc, &
                         c,r,delt,it,it0,dir,file12,    &
                         nroot,nrank,mnpr,opsum,ncomw,nerr)
     endif
  enddo

  call MPI_FINALIZE(nerr)


end program main
