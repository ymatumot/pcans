module fio

  implicit none

  private

  public :: fio__output
  public :: fio__input
  public :: fio__param
  public :: fio__energy

contains

  subroutine fio__output(up,uf,np,nx,nsp,np2,c,q,r,delt,delx,bc,it,it0,dir,file)

    integer, intent(in) :: np, nx, nsp, bc ,it, it0
    integer, intent(in) :: np2(1:nx+bc,nsp)
    real(8), intent(in) :: up(4,np,1:nx+bc,1,nsp)
    real(8), intent(in) :: uf(6,0:nx+1)
    real(8), intent(in) :: c,q(nsp), r(nsp), delt, delx
    character(len=*), intent(in) :: dir, file
    integer :: it2
    character(len=256) :: filename

    it2=it+it0

    !filename
    write(filename,'(a,i6.6,a)')trim(dir),it2,'_'//trim(file)
    open(10,file=filename,form='unformatted')

    !time & parameters
    write(10)it2,np,nx,nsp,bc,delt,delx,c
    write(10)np2
    write(10)q
    write(10)r

    !field data
    write(10)uf

    !particle data
    write(10)up

    
    close(10)

  end subroutine fio__output


  subroutine fio__input(up,uf,np2,c,q,r,delt,delx,it0,np,nx,nsp,bc,dir,file)

    integer, intent(in)  :: np, nsp, nx, bc
    character(len=*), intent(in) :: dir, file
    real(8), intent(out) :: up(4,np,1:nx+bc,nsp)
    real(8), intent(out) :: uf(6,0:nx+1)
    integer, intent(out) :: np2(1:nx+bc,nsp), it0
    real(8), intent(out) :: c, q(nsp), r(nsp), delt, delx
    integer :: inp, inx, insp, ibc

    !filename
    open(11,file=trim(dir)//trim(file),form='unformatted')

    !time & parameters
    read(11) it0,inp,inx,insp,ibc,delt,delx,c
    if((inx /= nx) .or. (inp /= np) .or. (insp /= nsp) .or. (ibc /= bc))then
       write(6,*) '** parameter mismatch **'
       stop
    endif
    read(11)np2
    read(11)q
    read(11)r

    !field data
    read(11)uf

    !particle data
    read(11)up

    close(11)

  end subroutine fio__input


  subroutine fio__param(np,nx,nsp,np2,c,q,r,vti,vte,va,rtemp,fpe,fge,ldb,delt,delx,bc,dir,file)

    integer, intent(in) :: np, nsp, nx, bc
    integer, intent(in) :: np2(1:nx+bc,nsp)
    real(8), intent(in) :: c, q(nsp), r(nsp), vti, vte, va, rtemp, fpe, fge, ldb, delt, delx
!    real(8), parameter  :: pi=4.0*atan(1.0)
    real(8), parameter  :: pi=3.141519
    character(len=*), intent(in) :: dir, file
    integer :: isp

    !filename
    open(9,file=trim(dir)//trim(file),status='unknown')

    write(9,610) nx, ldb
    write(9,620) (np2(1,isp),isp=1,nsp),np
    write(9,630) delx,delt,c
    write(9,640) (r(isp),isp=1,nsp)
    write(9,650) (q(isp),isp=1,nsp)
    write(9,660) fpe,fge,fpe*sqrt(r(2)/r(1)),fge*r(2)/r(1)
    write(9,670) va,vti,vte,(vti/va)**2,rtemp,vti/(fge*r(2)/r(1))
    write(9,*)
610 format(' grid size, debye lngth =======> ',i6,f8.4)
620 format(' particle number in cell=======> ',i6,i6,'/',i6)
630 format(' dx, dt, c ====================> ',f8.4,3x,f8.4,3x,f8.4)
640 format(' Mi, Me  ======================> ',2(1p,e10.2,1x))
650 format(' Qi, Qe  ======================> ',2(1p,e10.2,1x))
660 format(' Fpe, Fge, Fpi Fgi ============> ',4(1p,e10.2,1x))
670 format(' Va, Vi, Ve, beta, Te/Ti, rgi==> ',6(1p,e10.2,1x))
    close(9)
    
  end subroutine fio__param


  subroutine fio__energy(up,uf,np,nx,nsp,np2,c,r,delt,bc,it,it0,dir,file)

    integer, intent(in)          :: it, it0, np, nsp, nx, bc
    integer, intent(in)          :: np2(1:nx+bc,nsp)
    real(8), intent(in)          :: c, r(nsp), delt
    real(8), intent(in)          :: up(4,np,1:nx+bc,nsp)
    real(8), intent(in)          :: uf(6,0:nx+1)
    character(len=*), intent(in) :: dir, file
    integer :: i, ii, isp
    real(8) :: vene(nsp)
    real(8) :: efield, bfield, gam, total, u2
    real(8) :: pi

    pi = 4.0*atan(1.0)

    !filename
    if(it == 0) open(12,file=trim(dir)//trim(file),status='unknown')

    !energy
    vene(1:nsp) = 0.0
    efield = 0.0
    bfield = 0.0

    do isp=1,nsp
       do i=1,nx+bc
          do ii=1,np2(i,isp)
             u2 =  up(2,ii,i,isp)*up(2,ii,i,isp) &
                  +up(3,ii,i,isp)*up(3,ii,i,isp) &
                  +up(4,ii,i,isp)*up(4,ii,i,isp)
             gam = dsqrt(1.0+u2/(c*c))
             vene(isp) = vene(isp)+r(isp)*u2/(gam+1.)
          enddo
       enddo
    enddo
 
    do i=1,nx
       bfield = bfield+uf(1,i)*uf(1,i)+uf(2,i)*uf(2,i)+uf(3,i)*uf(3,i)
       efield = efield+uf(4,i)*uf(4,i)+uf(5,i)*uf(5,i)+uf(6,i)*uf(6,i)
    enddo
    efield = efield/(8.0*pi)
    bfield = bfield/(8.0*pi)

    total=vene(1)+vene(2)+efield+bfield

    write(12,610) (it+it0)*delt,vene(1),vene(2),efield,bfield,total
610 format(f8.2,5(e12.4))

  end subroutine fio__energy


end module fio
