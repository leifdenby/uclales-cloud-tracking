module modarray
  use tracking_common, only: cellptr

  private
  interface increase_array
    module procedure increase_array_i
    module procedure increase_array_r
    module procedure increase_array_p
  end interface increase_array

  public increase_array

  contains

  !> Increase an array of reals by allocating a new array of size `newsize` and copy the old array to the beginning of the new if
  !> the old one was allocated
  subroutine increase_array_r(array, newsize)
    real,    dimension(:,:), allocatable, intent(inout) :: array
    integer, dimension(:), intent(in)      :: newsize
    real,    dimension(:,:), allocatable   :: tmp

    if(allocated(array)) call move_alloc(array, tmp)
    allocate (array(newsize(1), newsize(2)))
    if (allocated(tmp)) then
      array(1:size(tmp,1),1:size(tmp,2)) = tmp
    end if

  end subroutine increase_array_r

  !> Increase an array of integers by allocating a new array of size `newsize` and copy the old array to the beginning of the new if
  !> the old one was allocated
  subroutine increase_array_i(array, newsize)
    integer, dimension(:,:),  allocatable, intent(inout) :: array
    integer, dimension(:), intent(in)      :: newsize
    integer, dimension(:,:), allocatable   :: tmp

    if(allocated(array)) call move_alloc(array, tmp)
    allocate (array(newsize(1), newsize(2)))
    if (allocated(tmp)) then
      array(1:size(tmp,1),1:size(tmp,2)) = tmp
    end if

  end subroutine increase_array_i

  !> Increase an array of the `cellptr` datatype by allocating a new array of size `newsize` and copy the old array to the beginning of the new if
  !> the old one was allocated
  subroutine increase_array_p(array, newsize)
    type(cellptr), dimension(:),  allocatable, intent(inout) :: array
    integer, intent(in)                        :: newsize
    type(cellptr), dimension(:),  allocatable  :: tmp

    if(allocated(array)) call move_alloc(array, tmp)
    allocate (array(newsize))
    if (allocated(tmp)) then
      array(1:size(tmp)) = tmp
    end if

  end subroutine increase_array_p
end module modarray
