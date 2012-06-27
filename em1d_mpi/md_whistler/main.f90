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
!    one-dimensional electromagnetic plasma simulation code
!
!    written by M Hoshino,  ISAS, 1984/09/12
!    revised  1985/03/08  1985/04/05  1997/05/06
!    revised for CANS (by Y. Matsumoto, STEL)  2004/06/22
!    re-written in F90 (by Y. Matsumoto, STEL)  2008/10/21
!    MPI parallelization (by Y. Matsumoto, STEL)  2009/4/1
!
!**********************************************************************c

  call init__set_param
  call fio__energy(up(2:4,1:np,nxs:nxe+bcp,1:nsp),uf,    &
                   np,nsp,np2,nxs,nxe,nxs1,nxe1,bcp, &
                   c,r,delt,it,it0,dir,file12,       &
                   nroot,nrank,mnpr,opsum,ncomw,nerr)
  call fio__output(up,uf,np,nxgs,nxge,nxs,nxe,nxs1,nxe1,nsp,np2,bc,bcp,nproc,nrank, &
                   c,q,r,delt,delx,0,it0,dir)

  do it=1,itmax-it0

     call particle__solv(gp,up,uf,                         &
                         c,q,r,delt,                       &
                         np,nsp,np2,nxs,nxe,nxs1,nxe1,bcp)
     call field__fdtd_i(uf,up,gp,                                      &
                        np,nsp,np2,nxgs,nxge,nxs,nxe,nxs1,nxe1,bc,bcp, &
                        q,c,delx,delt,gfac,                            &
                        nup,ndown,nroot,nproc,nrank,mnpr,opsum,nstat,ncomw,nerr)
     call boundary__particle(up,                                            &
                             np,nsp,np2,nxgs,nxge,nxs,nxe,nxs1,nxe1,bc,bcp, &
                             nup,ndown,nstat,mnpi,mnpr,ncomw,nerr)

     if(mod(it+it0,intvl1) == 0)                                                            &
          call fio__output(up,uf,np,nxgs,nxge,nxs,nxe,nxs1,nxe1,nsp,np2,bc,bcp,nproc,nrank, &
                           c,q,r,delt,delx,it,it0,dir)
     if(mod(it+it0,intvl2) == 0)                                                        &
          call fio__energy(up(2:4,1:np,nxs:nxe+bcp,1:nsp),uf,                           &
                           np,nsp,np2,nxs,nxe,nxs1,nxe1,bcp,c,r,delt,it,it0,dir,file12, &
                           nroot,nrank,mnpr,opsum,ncomw,nerr)
  enddo

  call MPI_FINALIZE(nerr)


end program main
