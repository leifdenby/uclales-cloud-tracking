module tracking_common
  implicit none

  private

  type :: cellptr
    type(celltype), pointer :: p
  end type cellptr

  type :: celltype
    integer :: id
    integer :: n_points
    integer(kind=2), allocatable, dimension(:,:)     :: loc
    integer(kind=2), allocatable, dimension (:,:)       :: value
    integer                                  :: nsplitters
    type(cellptr), allocatable, dimension(:) :: splitters
    integer                                  :: nparents
    type(cellptr), allocatable, dimension(:) :: parents
    integer                                  :: nchildren
    type(cellptr), allocatable, dimension(:) :: children
    integer                                  :: cloudsystemnr
    integer                                  :: cloudtype
!     type(cellptr), allocatable, dimension(:) :: siblings
    type(celltype), pointer :: previous
    type(celltype), pointer :: next
    type(celltype), pointer :: head
  end type celltype

  ! for indexing into %value array of celltype
  integer, parameter :: ibase = 1, itop = 2, ivalue = 3

  !> Marks all points which are outside objects in object mask
  integer, parameter :: OUTSIDE_OBJECTS = -1
  !> Marks all points which are inside objects, but not yet processed in object mask
  integer, parameter :: INSIDE_OBJECTS = 0
  !> For marking an object currently being worked on in object mask
  integer, parameter :: UNCLASSIFIED_IN_OBJECT = -2
  !> For marking an object which has already been processed in object mask
  integer, parameter :: PROCESSED_OBJECT = -3

  integer :: nx = -1
  integer :: ny = -1
  integer :: tstart = -1
  integer :: tend = -1
  real    :: dx, dy, dt
  logical, parameter :: ldebug = .false., lsiblings = .false.

  character(100) :: simulation_id

  public simulation_id

  ! both
  public ibase, itop, ivalue
  public celltype, cellptr
  public dt, dx, dy
  public nx, ny
  public tstart, tend

  public createcell, deletecell, firstcell, nextcell

  public print_cell_debug

  public UNCLASSIFIED_IN_OBJECT, OUTSIDE_OBJECTS, INSIDE_OBJECTS, PROCESSED_OBJECT

  public count_num_cells

  contains

  subroutine print_cell_debug(cell)
    type(celltype), pointer, intent(in) :: cell
    print *, ""
    print *, "cell id=", cell%id, "n_points=", cell%n_points
    print *, "nsplitters=", cell%nsplitters
    print *, "nparents=", cell%nparents
    print *, "nchildren=", cell%nchildren
    print *, "cloudsystemnr=", cell%cloudsystemnr
    print *, "cloudtype=", cell%cloudtype
  end subroutine print_cell_debug

  !> Given the input `cell` create a new `celltype` instance with the same
  !>`head` and `next`, with `previous` pointing to the current cell and set
  !>`next` of the current `cell` to the new `celltype` instance.
  !
  !! @ returns (through reassignment) a pointer to the new `celltype` instance
  subroutine createcell(cell)
    type(celltype), pointer, intent(inout) :: cell
    type(celltype), pointer       :: tmp

    if (associated(cell)) then
      allocate(tmp)
      tmp%previous => cell
      tmp%next     => cell%next
      tmp%head     => cell%head
      cell%next    => tmp
      cell => tmp
    else
      allocate(cell)
      cell%head => cell
      cell%previous => null()
      cell%next     => null()
    end if

    cell%n_points = 0
    cell%nchildren = 0
    cell%nsplitters= 0
    cell%nparents  = 0
    cell%cloudsystemnr = 0
    cell%cloudtype = -1
  end subroutine createcell

  !> Given `cell` remove all it and all associated allocated "splitters",
  !> "parents", "children"
  !
  !! The `next`, `previous` linked-list pointers are re-assigned so that a new
  !! contigous link-list is created, and a pointer into that list (either the
  !! "next" or "previous" element) is re-assigned to the `cell` argument
  subroutine deletecell(cell)
    type(celltype), pointer, intent(inout) :: cell
    type(celltype), pointer       :: tmp, tmp2
    integer :: ierr

    ierr = 0

    tmp => NULL()
    if (associated(cell%previous)) then
      tmp  => cell%previous
      tmp%next=> cell%next
    elseif (associated(cell%next)) then
      ! if there is no previous element we are deleting the head of the
      ! linked-list, and so need to update the head of every element to point to
      ! the new head
      tmp  => cell%next
      tmp%previous => null()
      tmp2 => tmp
      do
         if (ierr == -1) then
            exit
         else
            tmp2%head => tmp
            ierr = nextcell(tmp2)
         endif
      end do
    end if
    deallocate(cell%loc, cell%value, cell%splitters, cell%parents, cell%children, stat=ierr)
    deallocate(cell)
    if (associated(tmp)) then
      cell=>tmp
    else
      cell => null()
    end if
  end subroutine deletecell

  !> Through using `deletecell` de-allocate all cells (and associated "parents",
  !> etc) in the `next`, `previous` linked-list that the argument `cell` exists in
  subroutine delete_all(cell)
    type(celltype), pointer, intent(inout) :: cell
    integer :: iret
    iret = firstcell(cell)
    iret = 0
    do
    iret = iret + 1
      if (.not. associated(cell)) exit
      call deletecell(cell)
    end do
  end subroutine delete_all

  !> Assign the `cell` argument to the "next" `celltype` attribute of the supplied cell if the
  !> supplied cell is pointing to something
  !! @returns 0: `cell` is associated, -1: otherwise
  integer function nextcell(cell)
    type(celltype), pointer, intent(inout) :: cell

    nextcell = -1
    if (associated(cell%next)) then
      cell => cell%next
      nextcell = 0
    end if
  end function nextcell

  !> Assign the `cell` argument to the "head" `celltype` attribute of the supplied cell if the
  !supplied cell is pointing to something
  !@ return 0: `cell` is associated, -1: otherwise
  integer function firstcell(cell)
    type(celltype), pointer, intent(inout) :: cell
    firstcell = -1
    if(associated(cell)) then
      cell => cell%head
      firstcell = 0
    end if
  end function firstcell

  integer function count_num_cells(cell) result(num_cells)
    type(celltype), pointer, intent(in) :: cell

    type(celltype), pointer       :: tmp
    integer :: iret

    num_cells = 0
    tmp => cell

    iret = firstcell(tmp)
    do
      if (iret == -1) exit
      iret = nextcell(tmp)
      num_cells = num_cells + 1
    end do
  end function count_num_cells
end module tracking_common
