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

  integer :: nrel_max
  integer :: nx = -1
  integer :: ny = -1
  integer :: nt = -1
  integer :: tstart = -1
  real    :: dx, dy, dt
  integer, parameter :: minparentel = 100
  integer(kind=2)    :: cbstep
  logical, parameter :: ldebug = .false., lsiblings = .false.

  character(100) :: simulation_id

  public simulation_id

  ! both
  public ibase, itop, ivalue
  public celltype, cellptr
  public dt, dx, dy
  public minparentel
  public nt, nx, ny
  public nrel_max
  public tstart
  public cbstep

end module tracking_common
