module fio

  implicit none

  private

  public :: fio__init
  public :: fio__output
  public :: fio__input
  public :: fio__param
  public :: fio__energy
  public :: fio__mom
  public :: fio__psd
  public :: fio__progress_bar

  logical, save              :: is_init = .false.
  integer, save              :: np, nx, nsp, bc
  real(8), save              :: c, delx, delt, pi
  real(8), save, allocatable :: q(:), r(:)
  character(len=128), save   :: dir, dir_mom, dir_psd


contains


  subroutine fio__init(npin,nxin,nspin,bcin,qin,rin,cin,delxin,deltin,piin,dirin,dir_momin,dir_psdin)

    integer, intent(in) :: npin, nxin, nspin, bcin
    real(8), intent(in) :: piin
    real(8), intent(in) :: qin(nspin), rin(nspin), cin, delxin, deltin
    character(len=*), intent(in) :: dirin, dir_momin, dir_psdin

    np  = npin
    nx  = nxin
    nsp = nspin
    bc  = bcin
    allocate(q(nsp))
    allocate(r(nsp))
    q    = qin
    r    = rin
    c    = cin
    delx = delxin
    delt = deltin
    pi   = piin
    dir  = dirin
    dir_mom = dir_momin
    dir_psd = dir_psdin
    is_init = .true.

  end subroutine fio__init


  subroutine fio__output(up,uf,np2,it,file)

    integer, intent(in) :: it
    integer, intent(in) :: np2(1:nx+bc,nsp)
    real(8), intent(in) :: up(4,np,1:nx+bc,nsp)
    real(8), intent(in) :: uf(6,0:nx+1)
    character(len=*), intent(in) :: file
    character(len=256) :: filename

    if(.not.is_init)then
       write(6,*)'Initialize first by calling fio__init()'
       stop
    endif

    !filename
    write(filename,'(a,i6.6,a)')trim(dir),it,'_'//trim(file)
    open(10,file=filename,form='unformatted')

    !time & parameters
    write(10)it,np,nx,nsp,bc,delt,delx,c
    write(10)np2
    write(10)q
    write(10)r

    !field data
    write(10)uf

    !particle data
    write(10)up

    close(10)

  end subroutine fio__output


  subroutine fio__input(up,uf,np2,it0,file)

    character(len=*), intent(in) :: file
    real(8), intent(out) :: up(4,np,1:nx+bc,nsp)
    real(8), intent(out) :: uf(6,0:nx+1)
    integer, intent(out) :: np2(1:nx+bc,nsp), it0
    integer :: inp, inx, insp, ibc
    real(8) :: idelx, idelt, ic

    if(.not.is_init)then
       write(6,*)'Initialize first by calling fio__init()'
       stop
    endif

    !filename
    open(11,file=trim(dir)//trim(file),form='unformatted')

    !time & parameters
    read(11) it0,inp,inx,insp,ibc,delt,delx,c
    if((inx /= nx) .or. (inp /= np) .or. (insp /= nsp) .or. (ibc /= bc) &
      .or. (idelt /= delt) .or. (idelx /= delx) .or. ic /= c)then
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


  subroutine fio__param(np2,vti,vte,va,rtemp,fpe,fge,ldb,file)

    integer, intent(in) :: np2(1:nx+bc,nsp)
    real(8), intent(in) :: vti, vte, va, rtemp, fpe, fge, ldb
    character(len=*), intent(in) :: file
    integer :: isp

    if(.not.is_init)then
       write(6,*)'Initialize first by calling fio__init()'
       stop
    endif

    !filename
    open(9,file=trim(dir)//trim(file),status='replace')
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


  subroutine fio__energy(up,uf,np2,it,file)

    integer, intent(in)          :: it
    integer, intent(in)          :: np2(1:nx+bc,nsp)
    real(8), intent(in)          :: up(4,np,1:nx+bc,nsp)
    real(8), intent(in)          :: uf(6,0:nx+1)
    character(len=*), intent(in) :: file
    integer :: i, ii, isp
    real(8) :: vene(nsp)
    real(8) :: efield, bfield, gam, total, u2

    if(.not.is_init)then
       write(6,*)'Initialize first by calling fio__init()'
       stop
    endif

    !filename
    if(it == 0)then
       open(12,file=trim(dir)//trim(file),status='replace')
    else
       open(12,file=trim(dir)//trim(file),position='append')
    endif

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
             gam = sqrt(1.0D0+u2/(c*c))
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

    write(12,610) it*delt,vene(1),vene(2),efield,bfield,total
610 format(f8.2,5(e12.4))

  end subroutine fio__energy


  subroutine fio__mom(den,vel,temp,uf,it)

    integer, intent(in)    :: it
    real(8), intent(inout) :: den(0:nx+1,nsp), vel(0:nx+1,3,nsp), temp(0:nx+1,3,nsp), uf(6,0:nx+1)
    integer :: i
    real(8) :: tmp(0:nx+1,6)
    character(len=256) :: filename

    write(filename,'(a,i6.6,a)')trim(dir_mom),it,'_'//'den_i.dat'
    open(10,file=filename,status='replace')
    write(filename,'(a,i6.6,a)')trim(dir_mom),it,'_'//'den_e.dat'
    open(11,file=filename,status='replace')
    write(filename,'(a,i6.6,a)')trim(dir_mom),it,'_'//'Txx_i.dat'
    open(12,file=filename,status='replace')
    write(filename,'(a,i6.6,a)')trim(dir_mom),it,'_'//'Tyy_i.dat'
    open(13,file=filename,status='replace')
    write(filename,'(a,i6.6,a)')trim(dir_mom),it,'_'//'Tzz_i.dat'
    open(14,file=filename,status='replace')
    write(filename,'(a,i6.6,a)')trim(dir_mom),it,'_'//'Txx_e.dat'
    open(15,file=filename,status='replace')
    write(filename,'(a,i6.6,a)')trim(dir_mom),it,'_'//'Tyy_e.dat'
    open(16,file=filename,status='replace')
    write(filename,'(a,i6.6,a)')trim(dir_mom),it,'_'//'Tzz_e.dat'
    open(17,file=filename,status='replace')

    write(filename,'(a,i6.6,a)')trim(dir_mom),it,'_'//'vxi.dat'
    open(18,file=filename,status='replace')
    write(filename,'(a,i6.6,a)')trim(dir_mom),it,'_'//'vyi.dat'
    open(19,file=filename,status='replace')
    write(filename,'(a,i6.6,a)')trim(dir_mom),it,'_'//'vzi.dat'
    open(20,file=filename,status='replace')
    write(filename,'(a,i6.6,a)')trim(dir_mom),it,'_'//'vxe.dat'
    open(21,file=filename,status='replace')
    write(filename,'(a,i6.6,a)')trim(dir_mom),it,'_'//'vye.dat'
    open(22,file=filename,status='replace')
    write(filename,'(a,i6.6,a)')trim(dir_mom),it,'_'//'vze.dat'
    open(23,file=filename,status='replace')

    write(filename,'(a,i6.6,a)')trim(dir_mom),it,'_'//'bx.dat'
    open(24,file=filename,status='replace')
    write(filename,'(a,i6.6,a)')trim(dir_mom),it,'_'//'by.dat'
    open(25,file=filename,status='replace')
    write(filename,'(a,i6.6,a)')trim(dir_mom),it,'_'//'bz.dat'
    open(26,file=filename,status='replace')

    write(filename,'(a,i6.6,a)')trim(dir_mom),it,'_'//'ex.dat'
    open(27,file=filename,status='replace')
    write(filename,'(a,i6.6,a)')trim(dir_mom),it,'_'//'ey.dat'
    open(28,file=filename,status='replace')
    write(filename,'(a,i6.6,a)')trim(dir_mom),it,'_'//'ez.dat'
    open(29,file=filename,status='replace')


    !fields at x=i+1/2
    do i=1,nx+bc
       tmp(i,1) = uf(1,i)
       tmp(i,2) = 0.5*(+uf(2,i)+uf(2,i+1))
       tmp(i,3) = 0.5*(+uf(3,i)+uf(3,i+1))
       tmp(i,4) = 0.5*(+uf(4,i)+uf(4,i+1))
       tmp(i,5) = uf(5,i)
       tmp(i,6) = uf(6,i)
    enddo

    !propose vel and temperature
    vel(1:nx+bc,1,1:nsp) = vel(1:nx+bc,1,1:nsp)/den(1:nx+bc,1:nsp)
    vel(1:nx+bc,2,1:nsp) = vel(1:nx+bc,2,1:nsp)/den(1:nx+bc,1:nsp)
    vel(1:nx+bc,3,1:nsp) = vel(1:nx+bc,3,1:nsp)/den(1:nx+bc,1:nsp)
    temp(1:nx+bc,1,1:nsp) = +temp(1:nx+bc,1,1:nsp)/den(1:nx+bc,1:nsp) &
                            -vel(1:nx+bc,1,1:nsp)**2
    temp(1:nx+bc,2,1:nsp) = +temp(1:nx+bc,2,1:nsp)/den(1:nx+bc,1:nsp) &
                            -vel(1:nx+bc,2,1:nsp)**2
    temp(1:nx+bc,3,1:nsp) = +temp(1:nx+bc,3,1:nsp)/den(1:nx+bc,1:nsp) &
                            -vel(1:nx+bc,3,1:nsp)**2

    write(10,99)(den(i,1),i=1,nx+bc)
    write(11,99)(den(i,2),i=1,nx+bc)
    write(12,99)(temp(i,1,1),i=1,nx+bc)
    write(13,99)(temp(i,2,1),i=1,nx+bc)
    write(14,99)(temp(i,3,1),i=1,nx+bc)
    write(15,99)(temp(i,1,2),i=1,nx+bc)
    write(16,99)(temp(i,2,2),i=1,nx+bc)
    write(17,99)(temp(i,3,2),i=1,nx+bc)
    write(18,99)(vel(i,1,1),i=1,nx+bc)
    write(19,99)(vel(i,2,1),i=1,nx+bc)
    write(20,99)(vel(i,3,1),i=1,nx+bc)
    write(21,99)(vel(i,1,2),i=1,nx+bc)
    write(22,99)(vel(i,2,2),i=1,nx+bc)
    write(23,99)(vel(i,3,2),i=1,nx+bc)
    write(24,99)(tmp(i,1),i=1,nx+bc)
    write(25,99)(tmp(i,2),i=1,nx+bc)
    write(26,99)(tmp(i,3),i=1,nx+bc)
    write(27,99)(tmp(i,4),i=1,nx+bc)
    write(28,99)(tmp(i,5),i=1,nx+bc)
    write(29,99)(tmp(i,6),i=1,nx+bc)
99  format(100000E15.5)

    close(10)
    close(11)
    close(12)
    close(13)
    close(14)
    close(15)
    close(16)
    close(17)
    close(18)
    close(19)
    close(20)
    close(21)
    close(22)
    close(23)
    close(24)
    close(25)
    close(26)
    close(27)
    close(28)
    close(29)

  end subroutine fio__mom


  subroutine fio__psd(up,np2,it)

    integer, intent(in) :: np2(1:nx+bc,nsp)
    integer, intent(in) :: it
    real(8), intent(in) :: up(4,np,1:nx+bc,nsp)
    integer             :: i, ii, isp, ieq
    character(len=256)  :: filename

    write(filename,'(a,i6.6,a)')trim(dir_psd),it,'_'//'psd_i.dat'
    open(90,file=filename,status='replace')
    write(filename,'(a,i6.6,a)')trim(dir_psd),it,'_'//'psd_e.dat'
    open(91,file=filename,status='replace')

    isp=1
    do i=1,nx+bc
       do ii=1,np2(i,isp)
          write(90,'(4e20.5)')(up(ieq,ii,i,isp),ieq=1,4)
       enddo
    enddo

    isp=2
    do i=1,nx+bc
       do ii=1,np2(i,isp)
          write(91,'(4e20.5)')(up(ieq,ii,i,isp),ieq=1,4)
       enddo
    enddo

  end subroutine fio__psd


  subroutine fio__progress_bar(it,itmax)
    integer, intent(in) :: it, itmax
    integer :: k
    character(len=1) :: back
    character(len=128) :: bar

    write(bar,'(i3,1a,1x,1a,32a,32a)') &
         & 100*it/itmax, '%', '|', &
         & ('*', k=1,32*it/itmax), &
         & (' ', k=1,32-32*it/itmax), &
         & '|'

    back = char(8)
    write(6,'(128a)',advance='no') trim(bar)
    write(6,'(128a)',advance='no') (back,k=1,len(trim(bar)))

    call flush(6)

  end subroutine fio__progress_bar


end module fio
