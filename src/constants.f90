module constants
  implicit none
  private

  logical, parameter :: use_runtime_checks = .true.

  real, parameter :: real_maxval = real(huge(1_2))

  real, parameter :: trcpath_min = -1.
  real, parameter :: trcpath_max = 10000.
  real, parameter :: trcpath_thres = 100.

  real, parameter :: lwpmin   = -1.  !> min liquid water path value (in [kg/m^2]) used for scaling
  real, parameter :: lwpmax   = 10.  !> max liquid water path value (in [kg/m^2]) used for scaling
  real, parameter :: lwpthres = 0.01  !> liquid water path threshold (in [kg/m^2]) used for identifying clouds

  real, parameter :: coremin  = -5.
  real, parameter :: coremax  = 5.
  real, parameter :: corethres = 0.6  !> buoyancy threshold (in [K]) below which passive outflow is identified

  real, parameter :: rwpmin   = -1.
  real, parameter :: rwpmax   = 10.
  real, parameter :: rwpthres = 0.01

  real, parameter :: distmin = -4000.
  real, parameter :: distmax = 4000.

  integer, parameter :: nchunk = 100 !< number of time-steps to load in each "chunk"

  integer, parameter :: nmincells_cloud  = 1
  integer, parameter :: nmincells        = 4

  integer, parameter :: n_minparentel = 10  !< number of cells in a "parent" element necessary to make the parent-child connection

  ! parameters controlling growth of regions when doing splitting
  integer, parameter :: n_growth_steps_min = 5

  real, parameter    :: cbstep = 300.

  !> Distance allowed above say thermals to the cloud-base height [m]
  real, parameter    :: parent_array_allowed_height_offset = 200.


  !!! Derived variables below

  ! center and range values for parameters
  real, parameter :: thermzero  = 0.5*(trcpath_max + trcpath_min)
  real, parameter :: thermrange = (trcpath_max - trcpath_min)

  real, parameter :: lwpzero  = 0.5*(lwpmax + lwpmin)
  real, parameter :: lwprange = (lwpmax - lwpmin)

  real, parameter :: corezero  = 0.5*(coremax + coremin)
  real, parameter :: corerange = (coremax - coremin)

  real, parameter :: rwpzero  = 0.5*(rwpmax + rwpmin)
  real, parameter :: rwprange = (rwpmax - rwpmin)

  real, parameter :: distzero  = 0.5*(distmax + distmin)
  real, parameter :: distrange = (distmax - distmin)

  integer(kind=2), parameter :: i_lwpthres   = (lwpthres - lwpzero)/lwprange*real_maxval
  integer(kind=2), parameter :: i_corethres  = (corethres - corezero)/corerange*real_maxval
  integer(kind=2), parameter :: i_trcpath_thres = (trcpath_thres - thermzero)/thermrange*real_maxval
  integer(kind=2), parameter :: i_rwpthres   = (rwpthres - rwpzero)/rwprange*real_maxval

  integer(kind=2), parameter :: cbstep_as_int = cbstep/distrange*real_maxval

  public nchunk
  public distzero, distrange, distmin, distmax
  public corezero, corerange, coremin, coremax
  public lwpzero, lwprange, lwpmin, lwpmax
  public rwpzero, rwprange, rwpmin, rwpmax
  public thermzero, thermrange, trcpath_min, trcpath_max

  public i_corethres, i_lwpthres, i_trcpath_thres, i_rwpthres

  ! tracking
  public n_minparentel

  ! splitting
  public n_growth_steps_min

  !main
  public rwpthres, lwpthres, corethres, trcpath_thres

  public nmincells, nmincells_cloud

  public cbstep_as_int
  public real_maxval
  public parent_array_allowed_height_offset

  public use_runtime_checks

end module constants
