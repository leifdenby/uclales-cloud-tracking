module tracking_common
  implicit none

  private

  type :: cellptr
    type(celltype), pointer :: p
  end type cellptr

  type :: celltype
    integer :: id
    integer :: nelements
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
  integer, parameter :: MARKED_OBJECT = -2
  !> For marking an object which has already been processed in object mask
  integer, parameter :: PROCESSED_OBJECT = -3

  integer :: nrel_max
  integer :: nx = -1
  integer :: ny = -1
  integer :: nt = -1
  integer :: tstart = -1
  real    :: dx, dy, dt
  integer(kind=2)    :: cbstep
  logical, parameter :: ldebug = .false., lsiblings = .false.

  character(100) :: simulation_id

  public simulation_id

  ! both
  public ibase, itop, ivalue
  public celltype, cellptr
  public dt, dx, dy
  public nt, nx, ny
  public nrel_max
  public tstart
  public cbstep

  public createcell, deletecell, firstcell, nextcell

  public print_cell_debug

  public MARKED_OBJECT, OUTSIDE_OBJECTS, INSIDE_OBJECTS, PROCESSED_OBJECT

  contains

  subroutine print_cell_debug(cell)
    type(celltype), pointer, intent(in) :: cell
    print *, ""
    print *, "cell id=", cell%id, "nelements=", cell%nelements
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

    cell%nelements = 0
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
    type(celltype), pointer       :: tmp
    integer :: ierr
    tmp => NULL()
    if (associated(cell%previous)) then
      tmp  => cell%previous
      tmp%next=> cell%next
    elseif (associated(cell%next)) then
      tmp  => cell%next
      tmp%previous=> cell%previous
    end if
    deallocate(cell%loc, cell%value, cell%splitters, cell%parents, cell%children, stat=ierr)
    deallocate(cell)
    if (associated(tmp)) then
      cell=>tmp
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
end module tracking_common
