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


  integer :: nrel_max
  integer :: nx, ny, nt, tstart
  integer(kind=4), dimension(:,:,:), allocatable :: bool
  integer(kind=2), dimension(:,:,:,:), allocatable :: var
  integer(kind=2), dimension(:,:,:), allocatable :: var_min
  integer(kind=2), dimension(:,:,:), allocatable :: var_max
  integer, parameter :: ibase = 1, itop = 2
  integer :: ivalue = 3
  real    :: dx, dy, dt
  integer :: minparentel
  integer(kind=2)    :: cbstep
  logical, parameter :: ldebug = .false., lsiblings = .false.

  character(100) :: simulation_id

  public simulation_id

  ! both
  public bool, var
  public ibase, itop, ivalue
  public celltype, cellptr
  public dt, dx, dy
  public minparentel
  public nt, nx, ny
  public nrel_max
  public tstart
  public cbstep

  public var_min
  public var_max
end module tracking_common
