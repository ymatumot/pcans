module boundary

  implicit none

  private

  public  :: boundary__init
  public  :: boundary__dfield
  public  :: boundary__particle
  public  :: boundary__curre
  public  :: boundary__phi

  logical, save :: is_init = .false.
  integer, save :: np, nsp, nxgs, nxge, nygs, nyge, nxs, nxe, nys, nye, bc
  integer, save :: nup, ndown, mnpi, mnpr, ncomw
  integer       :: nerr
  integer, allocatable :: nstat(:)
  real(8), save :: delx, d_delx
  real(8), save :: u0x_bnd(2), u0y_bnd(2), u0z_bnd(2), delx_in


contains


  subroutine boundary__init(np_in,nsp_in,&
                            nxgs_in,nxge_in,nygs_in,nyge_in,nxs_in,nxe_in,nys_in,nye_in,bc_in, &
                            nup_in,ndown_in,mnpi_in,mnpr_in,ncomw_in,nerr_in,nstat_in, &
                            delx_in,u0x,u0y,u0z)

    integer, intent(in) :: np_in, nsp_in
    integer, intent(in) :: nxgs_in, nxge_in, nygs_in, nyge_in, nxs_in, nxe_in, nys_in, nye_in, bc_in
    integer, intent(in) :: nup_in, ndown_in, mnpi_in, mnpr_in, ncomw_in
    integer, intent(in) :: nerr_in, nstat_in(:)
    real(8), intent(in) :: delx_in
    real(8), optional, intent(in) :: u0x(2), u0y(2), u0z(2)

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
    bc    = bc_in
    nup   = nup_in
    ndown = ndown_in
    mnpi  = mnpi_in
    mnpr  = mnpr_in
    ncomw = ncomw_in
    nerr  = nerr_in

    allocate(nstat(size(nstat_in)))
    nstat = nstat_in

    delx = delx_in
    d_delx = 1./delx

    if(present(u0x))then
      u0x_bnd = u0x
    else
      u0x_bnd = (/0.0D0, 0.0D0/)
    endif
    if(present(u0y))then
      u0y_bnd = u0y
    else
      u0y_bnd = (/0.0D0, 0.0D0/)
    endif
    if(present(u0z))then
      u0z_bnd = u0z
    else
      u0z_bnd = (/0.0D0, 0.0D0/)
    endif

    is_init = .true.

  end subroutine boundary__init


  subroutine boundary__particle(up,np2)

    integer, intent(inout) :: np2(nys:nye,nsp)
    real(8), intent(inout) :: up(5,np,nys:nye,nsp)

    integer :: j, ii, iii, isp, ipos, jpos
    integer :: cnt(nys-1:nye+1,nsp), cnt2(nys:nye,nsp), cnt_tmp
    integer :: flag(np,nys:nye,nsp)
    real(8) :: bff_ptcl(np*5,nys-1:nye+1,nsp)

    if(.not.is_init)then
       write(6,*)'Initialize first by calling boundary__init()'
       stop
    endif

    cnt(nys-1:nye+1,1:nsp) = 0
    cnt2(nys:nye,1:nsp) = 0
    flag(1:np,nys:nye,1:nsp) = 0

    do isp=1,nsp
       do j=nys,nye
          do ii=1,np2(j,isp)

             ipos = floor(up(1,ii,j,isp)*d_delx)
             jpos = floor(up(2,ii,j,isp)*d_delx)

             if(bc==0)then
                if(ipos <= nxgs-1)then
                   up(1,ii,j,isp) = up(1,ii,j,isp)+(nxge-nxgs+1)*delx
                endif
                if(ipos >= nxge+1)then
                   up(1,ii,j,isp) = up(1,ii,j,isp)-(nxge-nxgs+1)*delx
                endif
             else if(bc==-1)then
                if(ipos <= nxgs-1)then
                   up(1,ii,j,isp) = 2.0*delx*nxgs-up(1,ii,j,isp)
                   up(3,ii,j,isp) = 2.0*u0x_bnd(1)-up(3,ii,j,isp)
                   up(4,ii,j,isp) = 2.0*u0y_bnd(1)-up(4,ii,j,isp)
                   up(5,ii,j,isp) = 2.0*u0z_bnd(1)-up(5,ii,j,isp)
                endif
                if(ipos >= nxge)then
                   up(1,ii,j,isp) = 2.0*delx*nxge-up(1,ii,j,isp)
                   up(3,ii,j,isp) = 2.0*u0x_bnd(2)-up(3,ii,j,isp)
                   up(4,ii,j,isp) = 2.0*u0y_bnd(2)-up(4,ii,j,isp)
                   up(5,ii,j,isp) = 2.0*u0z_bnd(2)-up(5,ii,j,isp)
                endif
             else
                write(*,*)'choose bc=0 (periodic) or bc=-1 (reflective)'
                stop
             endif

             if(jpos /= j)then
                if(jpos <= nygs-1)then
                   up(2,ii,j,isp) = up(2,ii,j,isp)+(nyge-nygs+1)*delx
                endif
                if(jpos >= nyge+1)then
                   up(2,ii,j,isp) = up(2,ii,j,isp)-(nyge-nygs+1)*delx
                endif
                bff_ptcl(1+5*cnt(jpos,isp),jpos,isp) = up(1,ii,j,isp)
                bff_ptcl(2+5*cnt(jpos,isp),jpos,isp) = up(2,ii,j,isp)
                bff_ptcl(3+5*cnt(jpos,isp),jpos,isp) = up(3,ii,j,isp)
                bff_ptcl(4+5*cnt(jpos,isp),jpos,isp) = up(4,ii,j,isp)
                bff_ptcl(5+5*cnt(jpos,isp),jpos,isp) = up(5,ii,j,isp)
                cnt(jpos,isp) = cnt(jpos,isp)+1
                cnt2(j,isp) = cnt2(j,isp)+1
                flag(cnt2(j,isp),j,isp) = ii
             endif
          enddo
       enddo
    enddo

    do isp=1,nsp
       !transfer to rank-1
       call MPI_SENDRECV(cnt(nys-1,isp),1,mnpi,ndown,100, &
                         cnt_tmp       ,1,mnpi,nup  ,100, &
                         ncomw,nstat,nerr)
       call MPI_SENDRECV(bff_ptcl(1               ,nys-1,isp),5*cnt(nys-1,isp),mnpr,ndown,101, &
                         bff_ptcl(5*cnt(nye,isp)+1,nye  ,isp),5*cnt_tmp       ,mnpr,nup  ,101, &
                         ncomw,nstat,nerr)
       cnt(nye,isp) = cnt(nye,isp)+cnt_tmp

       !transfer to rank+1
       call MPI_SENDRECV(cnt(nye+1,isp),1,mnpi,nup  ,200, &
                         cnt_tmp       ,1,mnpi,ndown,200, &
                         ncomw,nstat,nerr)
       call MPI_SENDRECV(bff_ptcl(1               ,nye+1,isp),5*cnt(nye+1,isp),mnpr,nup  ,201, &
                         bff_ptcl(5*cnt(nys,isp)+1,nys  ,isp),5*cnt_tmp       ,mnpr,ndown,201, &
                         ncomw,nstat,nerr)
       cnt(nys,isp) = cnt(nys,isp)+cnt_tmp
    enddo

    do isp=1,nsp
       do j=nys,nye
          iii=0
          cnt_tmp = cnt2(j,isp)
          loop1 :do ii=1,cnt2(j,isp)
             if(cnt(j,isp) == 0)then
                if(np2(j,isp) < flag(ii,j,isp)) exit loop1
                do while(np2(j,isp) == flag(cnt_tmp,j,isp))
                   np2(j,isp) = np2(j,isp)-1
                   if(np2(j,isp) < flag(ii,j,isp)) exit loop1
                   cnt_tmp = cnt_tmp-1
                enddo
                up(1:5,flag(ii,j,isp),j,isp) = up(1:5,np2(j,isp),j,isp)
                np2(j,isp) = np2(j,isp)-1
             else
                up(1,flag(ii,j,isp),j,isp) = bff_ptcl(1+5*iii,j,isp)
                up(2,flag(ii,j,isp),j,isp) = bff_ptcl(2+5*iii,j,isp)
                up(3,flag(ii,j,isp),j,isp) = bff_ptcl(3+5*iii,j,isp)
                up(4,flag(ii,j,isp),j,isp) = bff_ptcl(4+5*iii,j,isp)
                up(5,flag(ii,j,isp),j,isp) = bff_ptcl(5+5*iii,j,isp)
                iii = iii+1
                cnt(j,isp) = cnt(j,isp)-1
             endif
          enddo loop1
          
          if(cnt(j,isp) > 0)then
             do ii=1,cnt(j,isp)
                up(1,np2(j,isp)+ii,j,isp) = bff_ptcl(5*iii+1+5*(ii-1),j,isp)
                up(2,np2(j,isp)+ii,j,isp) = bff_ptcl(5*iii+2+5*(ii-1),j,isp)
                up(3,np2(j,isp)+ii,j,isp) = bff_ptcl(5*iii+3+5*(ii-1),j,isp)
                up(4,np2(j,isp)+ii,j,isp) = bff_ptcl(5*iii+4+5*(ii-1),j,isp)
                up(5,np2(j,isp)+ii,j,isp) = bff_ptcl(5*iii+5+5*(ii-1),j,isp)
             enddo
          endif
       enddo
    enddo

    do isp=1,nsp
       do j=nys,nye
          np2(j,isp) = np2(j,isp)+cnt(j,isp)
          if(np2(j,isp) > np) then
             write(*,*)"memory over (np2 > np)",np,np2(j,isp),j,isp
             stop
          endif
       enddo
    enddo

  end subroutine boundary__particle


  subroutine boundary__dfield(df)

    real(8), intent(inout) :: df(6,nxs-2:nxe+2,nys-2:nye+2)
    integer                :: i, ii
    real(8)                :: bff_snd(12*(nxe-nxs+1)), bff_rcv(12*(nxe-nxs+1))

    do i=nxs,nxe
       ii = 12*(i-nxs)
       bff_snd(ii+1)  = df(1,i,nys)
       bff_snd(ii+2)  = df(2,i,nys)
       bff_snd(ii+3)  = df(3,i,nys)
       bff_snd(ii+4)  = df(4,i,nys)
       bff_snd(ii+5)  = df(5,i,nys)
       bff_snd(ii+6)  = df(6,i,nys)
       bff_snd(ii+7)  = df(1,i,nys+1)
       bff_snd(ii+8)  = df(2,i,nys+1)
       bff_snd(ii+9)  = df(3,i,nys+1)
       bff_snd(ii+10) = df(4,i,nys+1)
       bff_snd(ii+11) = df(5,i,nys+1)
       bff_snd(ii+12) = df(6,i,nys+1)
    enddo
    call MPI_SENDRECV(bff_snd(1),12*(nxe-nxs+1),mnpr,ndown,110, &
                      bff_rcv(1),12*(nxe-nxs+1),mnpr,nup  ,110, &
                      ncomw,nstat,nerr)
    do i=nxs,nxe
       ii = 12*(i-nxs)
       df(1,i,nye+1) = bff_rcv(ii+1)   
       df(2,i,nye+1) = bff_rcv(ii+2)
       df(3,i,nye+1) = bff_rcv(ii+3)
       df(4,i,nye+1) = bff_rcv(ii+4)   
       df(5,i,nye+1) = bff_rcv(ii+5)
       df(6,i,nye+1) = bff_rcv(ii+6)
       df(1,i,nye+2) = bff_rcv(ii+7)   
       df(2,i,nye+2) = bff_rcv(ii+8)
       df(3,i,nye+2) = bff_rcv(ii+9)
       df(4,i,nye+2) = bff_rcv(ii+10)   
       df(5,i,nye+2) = bff_rcv(ii+11)
       df(6,i,nye+2) = bff_rcv(ii+12)
    enddo

    do i=nxs,nxe
       ii = 12*(i-nxs)
       bff_snd(ii+1)  = df(1,i,nye-1)
       bff_snd(ii+2)  = df(2,i,nye-1)
       bff_snd(ii+3)  = df(3,i,nye-1)
       bff_snd(ii+4)  = df(4,i,nye-1)
       bff_snd(ii+5)  = df(5,i,nye-1)
       bff_snd(ii+6)  = df(6,i,nye-1)
       bff_snd(ii+7)  = df(1,i,nye)
       bff_snd(ii+8)  = df(2,i,nye)
       bff_snd(ii+9)  = df(3,i,nye)
       bff_snd(ii+10) = df(4,i,nye)
       bff_snd(ii+11) = df(5,i,nye)
       bff_snd(ii+12) = df(6,i,nye)
    enddo
    call MPI_SENDRECV(bff_snd(1),12*(nxe-nxs+1),mnpr,nup  ,100, &
                      bff_rcv(1),12*(nxe-nxs+1),mnpr,ndown,100, &
                      ncomw,nstat,nerr)
    do i=nxs,nxe
       ii = 12*(i-nxs)
       df(1,i,nys-2) = bff_rcv(ii+1)   
       df(2,i,nys-2) = bff_rcv(ii+2)
       df(3,i,nys-2) = bff_rcv(ii+3)
       df(4,i,nys-2) = bff_rcv(ii+4)   
       df(5,i,nys-2) = bff_rcv(ii+5)
       df(6,i,nys-2) = bff_rcv(ii+6)
       df(1,i,nys-1) = bff_rcv(ii+7)   
       df(2,i,nys-1) = bff_rcv(ii+8)
       df(3,i,nys-1) = bff_rcv(ii+9)
       df(4,i,nys-1) = bff_rcv(ii+10)   
       df(5,i,nys-1) = bff_rcv(ii+11)
       df(6,i,nys-1) = bff_rcv(ii+12)
    enddo

    if(bc == 0)then
       df(1:6,nxs-2,nys-2:nye+2) = df(1:6,nxe-1,nys-2:nye+2)
       df(1:6,nxs-1,nys-2:nye+2) = df(1:6,nxe  ,nys-2:nye+2)
       df(1:6,nxe+1,nys-2:nye+2) = df(1:6,nxs  ,nys-2:nye+2)
       df(1:6,nxe+2,nys-2:nye+2) = df(1:6,nxs+1,nys-2:nye+2)
    else if(bc == -1)then
       df(1  ,nxs-2,nys-2:nye+2) = -df(1  ,nxs+1,nys-2:nye+2)
       df(2:4,nxs-2,nys-2:nye+2) = +df(2:4,nxs+2,nys-2:nye+2)
       df(5:6,nxs-2,nys-2:nye+2) = -df(5:6,nxs+1,nys-2:nye+2)

       df(1  ,nxs-1,nys-2:nye+2) = -df(1  ,nxs  ,nys-2:nye+2)
       df(2:4,nxs-1,nys-2:nye+2) = +df(2:4,nxs+1,nys-2:nye+2)
       df(5:6,nxs-1,nys-2:nye+2) = -df(5:6,nxs  ,nys-2:nye+2)

       df(1  ,nxe  ,nys-2:nye+2) = -df(1  ,nxe-1,nys-2:nye+2)
       df(2:4,nxe+1,nys-2:nye+2) = +df(2:4,nxe-1,nys-2:nye+2)
       df(5:6,nxe  ,nys-2:nye+2) = -df(5:6,nxe-1,nys-2:nye+2)

       df(1  ,nxe+1,nys-2:nye+2) = -df(1  ,nxe-2,nys-2:nye+2)
       df(2:4,nxe+2,nys-2:nye+2) = +df(2:4,nxe-2,nys-2:nye+2)
       df(5:6,nxe+1,nys-2:nye+2) = -df(5:6,nxe-2,nys-2:nye+2)
    else
       write(*,*)'choose bc=0 (periodic) or bc=-1 (reflective)'
       stop
    endif

  end subroutine boundary__dfield


  subroutine boundary__curre(uj)

    real(8), intent(inout) :: uj(3,nxs-2:nxe+2,nys-2:nye+2)
    integer                :: i, ii
    real(8)                :: bff_rcv(6*(nxe-nxs+4+1)), bff_snd(6*(nxe-nxs+4+1))

    do i=nxs-2,nxe+2
       ii = 6*(i-(nxs-2))
       bff_snd(ii+1) = uj(1,i,nys-2)
       bff_snd(ii+2) = uj(2,i,nys-2)
       bff_snd(ii+3) = uj(3,i,nys-2)
       bff_snd(ii+4) = uj(1,i,nys-1)
       bff_snd(ii+5) = uj(2,i,nys-1)
       bff_snd(ii+6) = uj(3,i,nys-1)
    enddo
    call MPI_SENDRECV(bff_snd(1),6*(nxe-nxs+4+1),mnpr,ndown,110, &
                      bff_rcv(1),6*(nxe-nxs+4+1),mnpr,nup  ,110, &
                      ncomw,nstat,nerr)
    do i=nxs-2,nxe+2
       ii = 6*(i-(nxs-2))
       uj(1,i,nye-1) = uj(1,i,nye-1)+bff_rcv(ii+1)
       uj(2,i,nye-1) = uj(2,i,nye-1)+bff_rcv(ii+2)
       uj(3,i,nye-1) = uj(3,i,nye-1)+bff_rcv(ii+3)
       uj(1,i,nye  ) = uj(1,i,nye  )+bff_rcv(ii+4)
       uj(2,i,nye  ) = uj(2,i,nye  )+bff_rcv(ii+5)
       uj(3,i,nye  ) = uj(3,i,nye  )+bff_rcv(ii+6)
    enddo

    do i=nxs-2,nxe+2
       ii = 6*(i-(nxs-2))
       bff_snd(ii+1) = uj(1,i,nye+1)
       bff_snd(ii+2) = uj(2,i,nye+1)
       bff_snd(ii+3) = uj(3,i,nye+1)
       bff_snd(ii+4) = uj(1,i,nye+2)
       bff_snd(ii+5) = uj(2,i,nye+2)
       bff_snd(ii+6) = uj(3,i,nye+2)
    enddo
    call MPI_SENDRECV(bff_snd(1),6*(nxe-nxs+4+1),mnpr,nup  ,100, &
                      bff_rcv(1),6*(nxe-nxs+4+1),mnpr,ndown,100, &
                      ncomw,nstat,nerr)
    do i=nxs-2,nxe+2
       ii = 6*(i-(nxs-2))
       uj(1,i,nys  ) = uj(1,i,nys  )+bff_rcv(ii+1)
       uj(2,i,nys  ) = uj(2,i,nys  )+bff_rcv(ii+2)
       uj(3,i,nys  ) = uj(3,i,nys  )+bff_rcv(ii+3)
       uj(1,i,nys+1) = uj(1,i,nys+1)+bff_rcv(ii+4)
       uj(2,i,nys+1) = uj(2,i,nys+1)+bff_rcv(ii+5)
       uj(3,i,nys+1) = uj(3,i,nys+1)+bff_rcv(ii+6)
    enddo

    do i=nxs-2,nxe+2
       ii = 3*(i-(nxs-2))
       bff_snd(ii+1) = uj(1,i,nys)
       bff_snd(ii+2) = uj(2,i,nys)
       bff_snd(ii+3) = uj(3,i,nys)
    enddo
    call MPI_SENDRECV(bff_snd(1),3*(nxe-nxs+4+1),mnpr,ndown,130, &
                      bff_rcv(1),3*(nxe-nxs+4+1),mnpr,nup  ,130, &
                      ncomw,nstat,nerr)
    do i=nxs-2,nxe+2
       ii = 3*(i-(nxs-2))
       uj(1,i,nye+1) = bff_rcv(ii+1)
       uj(2,i,nye+1) = bff_rcv(ii+2)
       uj(3,i,nye+1) = bff_rcv(ii+3)
    enddo

    do i=nxs-2,nxe+2
       ii = 3*(i-(nxs-2))
       bff_snd(ii+1) = uj(1,i,nye)
       bff_snd(ii+2) = uj(2,i,nye)
       bff_snd(ii+3) = uj(3,i,nye)
    enddo
    call MPI_SENDRECV(bff_snd(1),3*(nxe-nxs+4+1),mnpr,nup  ,140, &
                      bff_rcv(1),3*(nxe-nxs+4+1),mnpr,ndown,140, &
                      ncomw,nstat,nerr)
    do i=nxs-2,nxe+2
       ii = 3*(i-(nxs-2))
       uj(1,i,nys-1) = bff_rcv(ii+1)
       uj(2,i,nys-1) = bff_rcv(ii+2)
       uj(3,i,nys-1) = bff_rcv(ii+3)
    enddo

    if(bc == 0)then
       uj(1:3,nxs  ,nys-2:nye+2) = uj(1:3,nxs  ,nys-2:nye+2)+uj(1:3,nxe+1,nys-2:nye+2)
       uj(1:3,nxs+1,nys-2:nye+2) = uj(1:3,nxs+1,nys-2:nye+2)+uj(1:3,nxe+2,nys-2:nye+2)
       uj(1:3,nxe-1,nys-2:nye+2) = uj(1:3,nxe-1,nys-2:nye+2)+uj(1:3,nxs-2,nys-2:nye+2)
       uj(1:3,nxe  ,nys-2:nye+2) = uj(1:3,nxe  ,nys-2:nye+2)+uj(1:3,nxs-1,nys-2:nye+2)

       uj(1:3,nxs-1,nys-2:nye+2) = uj(1:3,nxe,nys-2:nye+2)
       uj(1:3,nxe+1,nys-2:nye+2) = uj(1:3,nxs,nys-2:nye+2)
    else if(bc == -1)then
       uj(1  ,nxs+1,nys-2:nye+2) = uj(1  ,nxs+1,nys-2:nye+2)-uj(1  ,nxs-1,nys-2:nye+2)
       uj(2:3,nxs  ,nys-2:nye+2) = uj(2:3,nxs  ,nys-2:nye+2)+uj(2:3,nxs-1,nys-2:nye+2)
       uj(2:3,nxs+1,nys-2:nye+2) = uj(2:3,nxs+1,nys-2:nye+2)+uj(2:3,nxs-2,nys-2:nye+2)
       uj(2:3,nxe-2,nys-2:nye+2) = uj(2:3,nxe-2,nys-2:nye+2)+uj(2:3,nxe+1,nys-2:nye+2)
       uj(2:3,nxe-1,nys-2:nye+2) = uj(2:3,nxe-1,nys-2:nye+2)+uj(2:3,nxe  ,nys-2:nye+2)
       uj(1  ,nxe-1,nys-2:nye+2) = uj(1  ,nxe-1,nys-2:nye+2)-uj(1  ,nxe+1,nys-2:nye+2)

       uj(1  ,nxs-1,nys-2:nye+2) = +uj(1  ,nxs+1,nys-2:nye+2)
       uj(2:3,nxs-1,nys-2:nye+2) = -uj(2:3,nxs  ,nys-2:nye+2)
       uj(2:3,nxe  ,nys-2:nye+2) = -uj(2:3,nxe-1,nys-2:nye+2)
       uj(1  ,nxe+1,nys-2:nye+2) = +uj(1  ,nxe-1,nys-2:nye+2)
    else
       write(*,*)'choose bc=0 (periodic) or bc=-1 (reflective)'
       stop
    endif

  end subroutine boundary__curre


  subroutine boundary__phi(phi,l)

    integer, intent(in)    :: l
    real(8), intent(inout) :: phi(nxs-1:nxe+1,nys-1:nye+1)
    integer                :: i, j, ii, bcc
    real(8)                :: bff_snd(nxe-nxs+1), bff_rcv(nxe-nxs+1)

    do i=nxs,nxe
       ii = i-nxs
       bff_snd(ii+1)  = phi(i,nys)
    enddo
    call MPI_SENDRECV(bff_snd(1),nxe-nxs+1,mnpr,ndown,110, &
                      bff_rcv(1),nxe-nxs+1,mnpr,nup  ,110, &
                      ncomw,nstat,nerr)
    do i=nxs,nxe
       ii = i-nxs
       phi(i,nye+1) = bff_rcv(ii+1)   
    enddo

    do i=nxs,nxe
       ii = i-nxs
       bff_snd(ii+1)  = phi(i,nye)
    enddo
    call MPI_SENDRECV(bff_snd(1),nxe-nxs+1,mnpr,nup  ,100, &
                      bff_rcv(1),nxe-nxs+1,mnpr,ndown,100, &
                      ncomw,nstat,nerr)
    do i=nxs,nxe
       ii = i-nxs
       phi(i,nys-1) = bff_rcv(ii+1)   
    enddo

    if(bc == 0)then
       phi(nxs-1,nys-1:nye+1) = phi(nxe,nys-1:nye+1)
       phi(nxe+1,nys-1:nye+1) = phi(nxs,nys-1:nye+1)
    else if(bc == -1)then
       select case(l)
          case(1)
          phi(nxs-1,nys-1:nye+1) = -phi(nxs  ,nys-1:nye+1)
          phi(nxe  ,nys-1:nye+1) = -phi(nxe-1,nys-1:nye+1)

          case(2,3)
          phi(nxs-1,nys-1:nye+1) = phi(nxs+1,nys-1:nye+1)
          phi(nxe+1,nys-1:nye+1) = phi(nxe-1,nys-1:nye+1)
       endselect
    else
       write(*,*)'choose bc=0 (periodic) or bc=-1 (reflective)'
       stop
    endif

  end subroutine boundary__phi


end module boundary
