module fio

  implicit none

  private

  public :: fio__init
  public :: fio__output
  public :: fio__input
  public :: fio__param
  public :: fio__energy


  logical, save :: is_init = .false.
  integer, save :: np, nsp, nxgs, nxge, nygs, nyge, nxs, nxe, nys, nye, nsfo, bc, n0
  integer, save :: nproc, nrank, nroot, nup, ndown, mnpr, opsum, ncomw
  integer       :: nerr
  integer, allocatable :: nstat(:)
  real(8), parameter :: pi = 4.0*atan(1.0d0)
  real(8), save :: c, delx, delt, gfac, d_delx, d_delt
  real(8), allocatable, save :: q(:), r(:)
  character(len=64), save :: dir, fname_param, fname_energy


contains


  subroutine fio__init(np_in,nsp_in,&
                       nxgs_in,nxge_in,nygs_in,nyge_in,nxs_in,nxe_in,nys_in,nye_in,nsfo_in,bc_in, &
                       q_in,r_in,c_in,delx_in,delt_in,gfac_in,n0_in,&
                       fname_param_in,fname_energy_in,dir_in, &
                       nproc_in,nrank_in,nroot_in,nup_in,ndown_in,mnpr_in,opsum_in,ncomw_in,nerr_in,nstat_in)

    integer, intent(in) :: np_in, nsp_in
    integer, intent(in) :: nxgs_in, nxge_in, nygs_in, nyge_in, nxs_in, nxe_in, nys_in, nye_in, nsfo_in, bc_in
    integer, intent(in) :: nproc_in, nroot_in, nrank_in, nup_in, ndown_in, mnpr_in, opsum_in, ncomw_in
    integer, intent(in) :: nerr_in, nstat_in(:)
    integer, intent(in) :: n0_in
    real(8), intent(in) :: q_in(nsp_in), r_in(nsp_in), c_in, delx_in, delt_in, gfac_in
    character(len=*), intent(in) :: dir_in, fname_energy_in, fname_param_in

    np    = np_in
    nsp   = nsp_in
    nxgs  = nxgs_in
    nxge  = nxge_in
    nygs  = nygs_in
    nyge  = nyge_in
    nxs   = nxs_in
    nxe   = nxe_in
    nys   = nys_in
    nye   = nye_in
    nsfo  = nsfo_in
    bc    = bc_in
    nproc = nproc_in
    nrank = nrank_in
    nroot = nroot_in
    nup   = nup_in
    ndown = ndown_in
    mnpr  = mnpr_in
    opsum = opsum_in
    ncomw = ncomw_in
    nerr  = nerr_in
    c     = c_in
    delx  = delx_in
    delt  = delt_in
    gfac  = gfac_in
    d_delx = 1./delx
    d_delt = 1./delt
    n0    = n0_in
    dir   = dir_in
    fname_param  = fname_param_in
    fname_energy = fname_energy_in

    allocate(nstat(size(nstat_in)))
    nstat = nstat_in
    allocate(q(nsp))
    q     = q_in
    allocate(r(nsp))
    r     = r_in

    is_init = .true.

  end subroutine fio__init


  subroutine fio__output(up,uf,np2,it)

    integer, intent(in) :: np2(nys:nye,nsp)
    integer, intent(in) :: it
    real(8), intent(in) :: up(5,np,nys:nye,nsp)
    real(8), intent(in) :: uf(6,nxs-2:nxe+2,nys-2:nye+2)
    character(len=256) :: filename

    if(.not.is_init)then
       write(6,*)'Initialize first by calling fio__init()'
       stop
    endif

    !filename
    write(filename,'(a,i6.6,a,i3.3,a)')trim(dir),it,'_rank=',nrank,'.dat'
    open(110+nrank,file=filename,form='unformatted')

    !time & parameters
    write(110+nrank)it,np,nxgs,nxge,nygs,nyge,nxs,nxe,nys,nye,nsp,nproc,bc,delt,delx,c
    write(110+nrank)np2
    write(110+nrank)q
    write(110+nrank)r

    !field data
    write(110+nrank)uf(1:6,nxs-1:nxe+1,nys-1:nye+1)

    !particle data
    write(110+nrank)up
       
    close(110+nrank)

  end subroutine fio__output


  subroutine fio__input(up,uf,np2,it0,file)

    character(len=*), intent(in) :: file
    integer, intent(out) :: np2(nys:nye,nsp), it0
    real(8), intent(out) :: up(5,np,nys:nye,nsp)
    real(8), intent(out) :: uf(6,nxs-2:nxe+2,nys-2:nye+2)
    integer :: inp, inxgs, inxge, inygs, inyge, inxs, inxe, inys, inye, insp, ibc, inproc
    real(8) :: ic, idelt, idelx

    if(.not.is_init)then
       write(6,*)'Initialize first by calling fio__init()'
       stop
    endif

    !filename
    open(110+nrank,file=trim(dir)//trim(file),form='unformatted')

    !time & parameters
    read(110+nrank)it0,inp,inxgs,inxge,inygs,inyge,inxs,inxe,inys,inye,insp,inproc,ibc,idelt,idelx,ic
    if((inxgs /= nxgs) .or. (inxge /= nxge)  .or.(inygs /= nygs) .or. (inyge /= nyge) &
       .or. (inxs /= nxs) .or. (inxe /= nxe)  .or.(inys /= nys) .or. (inye /= nye)    &
       .or. (inp /= np) .or. (insp /= nsp) .or. (ibc /= bc) .or. (idelt /= delt) .or. (idelx /= delx) &
       .or. (ic /= c) .or. (inproc /= nproc))then
       write(6,*) '** parameter mismatch **'
       stop
    endif
    read(110+nrank)np2
    read(110+nrank)q
    read(110+nrank)r

    !field data
    read(110+nrank)uf(1:6,nxs-1:nxe+1,nys-1:nye+1)

    !particle data
    read(110+nrank)up

    close(110+nrank)

  end subroutine fio__input


  subroutine fio__param(np2,temp,rtemp,fpe,fge,ldb)

    integer, intent(in) :: np2(nys:nye,nsp)
    real(8), intent(in) :: temp, rtemp, fpe, fge, ldb
    integer :: isp
    real(8) :: vti, vte, va

    if(.not.is_init)then
       write(6,*)'Initialize first by calling fio__init()'
       stop
    endif

    vti = sqrt(2.*temp/r(1))
    vte = sqrt(2.*temp*rtemp/r(2))
    va  = fge*r(2)*c/q(1)/sqrt(4.*pi*r(1)*n0)

    if(nrank == nroot)then

       !filename
       open(9,file=trim(dir)//trim(fname_param),status='replace')

       write(9,610) nxge-nxgs+1,'x',nyge-nygs+1,ldb
       write(9,620) (np2(nys,isp),isp=1,nsp),np
       write(9,630) delx,delt,c
       write(9,640) (r(isp),isp=1,nsp)
       write(9,650) (q(isp),isp=1,nsp)
       write(9,660) fpe,fge,fpe*sqrt(r(2)/r(1)),fge*r(2)/r(1)
       write(9,670) va,vti,vte,(vti/va)**2,rtemp,vti/(fge*r(2)/r(1))
       write(9,*)
610    format(' grid size, debye lngth =======> ',i6,a,i6,f8.4)
620    format(' particle number in cell ======> ',i8,i8,'/',i8)
630    format(' dx, dt, c ====================> ',f8.4,3x,f8.4,3x,f8.4)
640    format(' Mi, Me  ======================> ',2(1p,e10.2,1x))
650    format(' Qi, Qe  ======================> ',2(1p,e10.2,1x))
660    format(' Fpe, Fge, Fpi Fgi ============> ',4(1p,e10.2,1x))
670    format(' Va, Vi, Ve, beta, Te/Ti, rgi==> ',6(1p,e10.2,1x))
       close(9)

    endif
    
  end subroutine fio__param


  subroutine fio__energy(up,uf,np2,it)

    integer, intent(in) :: it,np2(nys:nye,nsp)
    real(8), intent(in) :: up(5,np,nys:nye,nsp)
    real(8), intent(in) :: uf(6,nxs-2:nxe+2,nys-2:nye+2)
    integer :: i, j, ii, isp
    real(8) :: vene(nsp), vene_g(nsp)
    real(8) :: efield, bfield, gam, total, u2
    real(8) :: efield_g, bfield_g

    if(.not.is_init)then
       write(6,*)'Initialize first by calling fio__init()'
       stop
    endif

    !filename
    if(nrank == nroot)then
       if(it == 0)then
          open(12,file=trim(dir)//trim(fname_energy),status='replace')
       else
          open(12,file=trim(dir)//trim(fname_energy),position='append')
       endif
    endif

    !energy
    vene(1:nsp) = 0.0
    efield = 0.0
    bfield = 0.0

    do isp=1,nsp
       do j=nys,nye
          do ii=1,np2(j,isp)
             u2 =  up(3,ii,j,isp)*up(3,ii,j,isp) &
                  +up(4,ii,j,isp)*up(4,ii,j,isp) &
                  +up(5,ii,j,isp)*up(5,ii,j,isp)
             gam = sqrt(1.0D0+u2/(c*c))
             vene(isp) = vene(isp)+r(isp)*u2/(gam+1.)
          enddo
       enddo
    enddo

    do isp=1,nsp
       call MPI_REDUCE(vene(isp),vene_g(isp),1,mnpr,opsum,nroot,ncomw,nerr)
    enddo

    do j=nys,nye
    do i=nxs,nxe+bc
       bfield = bfield+uf(1,i,j)*uf(1,i,j)+uf(2,i,j)*uf(2,i,j)+uf(3,i,j)*uf(3,i,j)
       efield = efield+uf(4,i,j)*uf(4,i,j)+uf(5,i,j)*uf(5,i,j)+uf(6,i,j)*uf(6,i,j)
    enddo
    enddo
    if(bc == -1)then
       i=nxe
       do j=nys,nye
          bfield = bfield+uf(2,i,j)*uf(2,i,j)+uf(3,i,j)*uf(3,i,j)
          efield = efield+uf(4,i,j)*uf(4,i,j)
       enddo
    endif

    efield = efield/(8.0*pi)
    bfield = bfield/(8.0*pi)
    call MPI_REDUCE(efield,efield_g,1,mnpr,opsum,nroot,ncomw,nerr)
    call MPI_REDUCE(bfield,bfield_g,1,mnpr,opsum,nroot,ncomw,nerr)

    if(nrank == nroot)then
       total=vene_g(1)+vene_g(2)+efield_g+bfield_g
       
       write(12,610) it*delt,vene_g(1),vene_g(2),efield_g,bfield_g,total
610    format(f8.2,5(1p,e12.4))
    endif

  end subroutine fio__energy


end module fio
