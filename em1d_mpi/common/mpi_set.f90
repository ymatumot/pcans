module mpi_set

  implicit none
  private
  public :: mpi_set__init, MPI_WTIME

  include 'mpif.h'

  integer, public              :: nerr, ncomw, nsize, nrank
  integer, parameter, public   :: mnpi  = MPI_INTEGER
  integer, parameter, public   :: mnpr  = MPI_DOUBLE_PRECISION
  integer, parameter, public   :: mnpc  = MPI_CHARACTER
  integer, parameter, public   :: opsum = MPI_SUM
  integer, public              :: nstat(MPI_STATUS_SIZE)
  integer, public              :: nxs, nxs1, nxe, nxe1, nup, ndown

contains

  subroutine mpi_set__init(nxgs,nxge,bc,nproc)

    integer, intent(in) :: nxgs, nxge, bc, nproc
    integer             :: iwork1, iwork2

    !*********** Initialization for MPI  ***************!
    call MPI_INIT(nerr)
    ncomw = MPI_COMM_WORLD
    call MPI_COMM_SIZE(ncomw, nsize, nerr)
    call MPI_COMM_RANK(ncomw, nrank, nerr)

    if(nsize /= nproc) then
       call MPI_ABORT(ncomw, 9, nerr)
       call MPI_FINALIZE(nerr)
       stop '** proc number mismatch **'
    endif

    !start and end of loop counter
    iwork1 = (nxge-nxgs+1)/nsize
    iwork2 = mod(nxge-nxgs+1,nsize)
    nxs = nrank*iwork1+nxgs+min(nrank,iwork2)
    nxe = nxs+iwork1-1
    if(iwork2 > nrank) nxe = nxe+1

    !For MPI_SENDRECV
    nxs1  = nxs-1
    nxe1  = nxe+1
    nup   = nrank+1
    ndown = nrank-1

    if(bc == 0)then
       ! periodic boundary condition
       if(nrank == nsize-1) nup = 0 
       if(nrank == 0) ndown = nsize-1
    else if(bc == -1)then
       ! reflective boundary condition
       if(nrank == nsize-1) nup = MPI_PROC_NULL
       if(nrank == 0)     ndown = MPI_PROC_NULL
    else
       write(*,*)'choose bc=0 (periodic) or bc=-1 (reflective)'
       stop
    endif

  end subroutine mpi_set__init


end module mpi_set

