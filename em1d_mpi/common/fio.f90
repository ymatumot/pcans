module fio

  implicit none

  private

  public :: fio__output
  public :: fio__input
  public :: fio__param
  public :: fio__energy


contains


  subroutine fio__output(up,uf,np,nxgs,nxge,nxs,nxe,nxs1,nxe1,nsp,np2,bc,bcp,nproc,nrank, &
                         c,q,r,delt,delx,it,it0,dir)

    integer, intent(in) :: np, nxgs, nxge, nxs, nxe, nxs1, nxe1, nsp, bc, bcp
    integer, intent(in) :: np2(nxs:nxe+bcp,nsp)
    integer, intent(in) :: nproc, nrank
    integer, intent(in) :: it, it0
    real(8), intent(in) :: up(4,np,nxs:nxe+bcp,nsp)
    real(8), intent(in) :: uf(6,nxs1:nxe1)
    real(8), intent(in) :: c, q(nsp), r(nsp), delt, delx
    character(len=*), intent(in) :: dir
    integer :: it2
    character(len=256) :: filename

    it2=it+it0

    !filename
    write(filename,'(a,i6.6,a,i3.3,a)')trim(dir),it2,'_rank=',nrank,'.dat'
    open(100+nrank,file=filename,form='unformatted')

    !time & parameters
    write(100+nrank)it2,np,nxgs,nxge,nsp,nproc,bc,bcp,delt,delx,c
    write(100+nrank)np2
    write(100+nrank)q
    write(100+nrank)r

    !field data
    write(100+nrank)uf

    !particle data
    write(100+nrank)up
       
    close(100+nrank)

  end subroutine fio__output


  subroutine fio__input(up,uf,np2,c,q,r,delt,delx,it0, &
                        np,nxgs,nxge,nxs,nxe,nxs1,nxe1,nsp,bc,bcp,nproc,nrank,dir,file)

    integer, intent(in)  :: np, nxgs, nxge, nxs, nxe, nxs1, nxe1, nsp, bc, bcp, nproc, nrank
    character(len=*), intent(in) :: dir, file
    integer, intent(out) :: np2(nxs:nxe+bcp,nsp), it0
    real(8), intent(out) :: up(4,np,nxs:nxe+bcp,nsp)
    real(8), intent(out) :: uf(6,nxs1:nxe1)
    real(8), intent(out) :: c, q(nsp), r(nsp), delt, delx
    integer :: inp, inxgs, inxge, insp, inproc, ibc, ibcp

    !filename
    open(101+nrank,file=trim(dir)//trim(file),form='unformatted')

    !time & parameters
    read(101+nrank)it0,inp,inxgs,inxge,insp,inproc,ibc,ibcp,delt,delx,c
    if((inxgs /= nxgs) .or. (inxge /= nxge)  .or. &
       (inp /= np) .or. (insp /= nsp) .or. (inproc /= nproc) .or. &
       (ibc /= bc) .or. (ibcp /= bcp))then
       write(6,*) '** parameter mismatch **'
       stop
    endif
    read(101+nrank)np2
    read(101+nrank)q
    read(101+nrank)r

    !field data
    read(101+nrank)uf

    !particle data
    read(101+nrank)up

    close(101+nrank)

  end subroutine fio__input


  subroutine fio__param(np,nsp,np2,nxgs,nxge,nxs,nxe,nxs1,nxe1,bcp, &
                        c,q,r,vti,vte,va,rtemp,fpe,fge,             &
                        ldb,delt,delx,dir,file,                     &
                        nroot,nrank)

    integer, intent(in)          :: np, nsp 
    integer, intent(in)          :: nxgs, nxge, nxs, nxe, nxs1, nxe1, bcp, nroot, nrank
    integer, intent(in)          :: np2(nxs:nxe+bcp,nsp)
    real(8), intent(in)          :: c, q(nsp), r(nsp), vti, vte, va, rtemp, fpe, fge, ldb, delt, delx
    character(len=*), intent(in) :: dir, file
    integer :: isp
    real(8) :: pi

!    pi = 4.0*atan(1.0)
    pi = 3.141519

    if(nrank == nroot)then

       !filename
       open(9,file=trim(dir)//trim(file),status='unknown')

       write(9,610) nxge-nxgs+1, ldb
       write(9,620) (np2(nxgs,isp),isp=1,nsp),np
       write(9,630) delx,delt,c
       write(9,640) (r(isp),isp=1,nsp)
       write(9,650) (q(isp),isp=1,nsp)
       write(9,660) fpe,fge,fpe*sqrt(r(2)/r(1)),fge*r(2)/r(1)
       write(9,670) va,vti,vte,(vti/va)**2,rtemp,vti/(fge*r(2)/r(1))
       write(9,*)
610    format(' grid size, debye lngth =======> ',i6,f8.4)
620    format(' particle number in cell=======> ',i6,i6,'/',i6)
630    format(' dx, dt, c ====================> ',f8.4,3x,f8.4,3x,f8.4)
640    format(' Mi, Me  ======================> ',2(1p,e10.2,1x))
650    format(' Qi, Qe  ======================> ',2(1p,e10.2,1x))
660    format(' Fpe, Fge, Fpi Fgi ============> ',4(1p,e10.2,1x))
670    format(' Va, Vi, Ve, beta, Te/Ti, rgi==> ',6(1p,e10.2,1x))
       close(9)

    endif
    
  end subroutine fio__param


  subroutine fio__energy(uv,uf,np,nsp,np2,nxs,nxe,nxs1,nxe1,bcp, &
                         c,r,delt,it,it0,dir,file,               &
                         nroot,nrank,mnpr,opsum,ncomw,nerr)

    integer, intent(in)          :: nxs, nxe, nxs1, nxe1, bcp, nroot, nrank, mnpr, opsum, ncomw
    integer, intent(in)          :: it, it0, np, nsp, np2(nxs:nxe+bcp,nsp)
    integer, intent(inout)       :: nerr
    real(8), intent(in)          :: c, r(nsp), delt
    real(8), intent(in)          :: uv(3,np,nxs:nxe+bcp,nsp)
    real(8), intent(in)          :: uf(6,nxs1:nxe1)
    character(len=*), intent(in) :: dir, file
    integer :: i, ii, isp
    real(8) :: pi
    real(8) :: vene(nsp), vene_g(nsp)
    real(8) :: efield, bfield, gam, total, u2
    real(8) :: efield_g, bfield_g

!    pi = 4.0*atan(1.0)
    pi = 3.141519

    !filename
    if(nrank == nroot)then
       if(it == 0) open(12,file=trim(dir)//trim(file),status='unknown')
    endif

    !energy
    vene(1:nsp) = 0.0
    efield = 0.0
    bfield = 0.0

    do isp=1,nsp
       do i=nxs,nxe+bcp
          do ii=1,np2(i,isp)
             u2 =  uv(1,ii,i,isp)*uv(1,ii,i,isp) &
                  +uv(2,ii,i,isp)*uv(2,ii,i,isp) &
                  +uv(3,ii,i,isp)*uv(3,ii,i,isp)
             gam = dsqrt(1.0+u2/(c*c))
             vene(isp) = vene(isp)+r(isp)*u2/(gam+1.)
          enddo
       enddo
    enddo

    do isp=1,nsp
       call MPI_REDUCE(vene(isp),vene_g(isp),1,mnpr,opsum,nroot,ncomw,nerr)
    enddo

    do i=nxs,nxe
       bfield = bfield+uf(1,i)**2+uf(2,i)**2+uf(3,i)**2
       efield = efield+uf(4,i)**2+uf(5,i)**2+uf(6,i)**2
    enddo
    efield = efield/(8.0*pi)
    bfield = bfield/(8.0*pi)
    call MPI_REDUCE(efield,efield_g,1,mnpr,opsum,nroot,ncomw,nerr)
    call MPI_REDUCE(bfield,bfield_g,1,mnpr,opsum,nroot,ncomw,nerr)

    if(nrank == nroot)then
       total=vene_g(1)+vene_g(2)+efield_g+bfield_g
       
       write(12,610) (it+it0)*delt,vene_g(1),vene_g(2),efield_g,bfield_g,total
610    format(f8.2,5(e12.4))
    endif

  end subroutine fio__energy


end module fio
