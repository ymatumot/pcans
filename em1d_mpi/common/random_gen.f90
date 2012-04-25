module random_gen

  implicit none

  private

  public :: random_gen__init, random_gen__bm

  real(8), save :: pi


contains


  subroutine random_gen__init(nrank)

    integer, intent(in)  :: nrank
    integer              :: n
    integer, allocatable :: seed(:)

    call random_seed()
    call random_seed(size=n)
    allocate(seed(n))
    call random_seed(get=seed)
    seed(1:n) = seed(1:n)+nrank
    call random_seed(put=seed)
    deallocate(seed)

    pi = 4.*atan(1.D0)

  end subroutine random_gen__init


  subroutine random_gen__bm(r1,r2)
    !Box-Muller transform

    real(8), intent(out) :: r1, r2
    real(8)              :: aa, bb

    call random_number(aa)
    call random_number(bb)

    r1 = dsqrt(-2.*dlog(aa))*cos(2.*pi*bb)
    r2 = dsqrt(-2.*dlog(aa))*sin(2.*pi*bb)

  end subroutine random_gen__bm


end module random_gen

