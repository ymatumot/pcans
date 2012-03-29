subroutine random__init(ranflag,id)
  implicit none
  integer, intent(in)  :: ranflag
  integer, intent(in)  :: id
  integer :: nseed
  integer,allocatable  :: seed(:)
  integer :: numach

  select case(ranflag)
  case(1) 
     call random_seed()
     call random_seed(size=nseed)
     allocate(seed(nseed))
     call random_seed(get=seed)
     seed(1:nseed) = seed(1:nseed) + id
     call random_seed(put=seed)
     deallocate(seed)
!   case(2) 
!   case(3) 
  case default
  end select
  return
end subroutine random__init

subroutine random__generate(N,ran,ranflag)
  implicit none
  
  integer,intent(in)                 :: N
  integer,intent(in)                 :: ranflag
  real(8),intent(inout),dimension(N) :: ran

  integer :: i
  integer :: idum
  idum=3145239

  select case(ranflag)
  case(1)
     call random_number(ran)
!   case(2) 
!   case(3) 
  case default
     call random_number(ran)
  end select
  return
end subroutine Random__generate
        
        

