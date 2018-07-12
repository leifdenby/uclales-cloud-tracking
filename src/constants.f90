module constants
  implicit none
  private

  real, parameter :: real_maxval = real(huge(1_2))

  real, parameter :: thermmin = -1.
  real, parameter :: thermmax = 10000.
  real, parameter :: thermthres = 100.

  real, parameter :: lwpmin   = -1.
  real, parameter :: lwpmax   = 10.
  real, parameter :: lwpthres = 0.01

  real, parameter :: coremin  = -5.
  real, parameter :: coremax  = 5.
  real, parameter :: corethres = 0.6  !> buoyancy threshold (in [K]) below which passive outflow is identified

  real, parameter :: rwpmin   = -1.
  real, parameter :: rwpmax   = 10.
  real, parameter :: rwpthres = 0.01

  real, parameter :: distmin = -2000.
  real, parameter :: distmax = 2000.

  integer, parameter :: nchunk = 100 !< number of time-steps to load in each "chunk"

  integer, parameter :: nmincells_cloud  = 1
  integer, parameter :: nmincells        = 4

  integer, parameter :: n_minparentel = 10  !< number of cells in a "parent" element necessary to make the parent-child connection

  ! parameters controlling growth of regions when doing splitting
  integer, parameter :: n_growth_steps_min = 5

  real, parameter    :: cbstep = 300.


  !!! Derived variables below

  ! center and range values for parameters
  real, parameter :: thermzero  = 0.5*(thermmax + thermmin)
  real, parameter :: thermrange = (thermmax - thermmin)

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
  integer(kind=2), parameter :: i_thermthres = (thermthres - thermzero)/thermrange*real_maxval
  integer(kind=2), parameter :: i_rwpthres   = (rwpthres - rwpzero)/rwprange*real_maxval

  integer(kind=2), parameter :: cbstep_as_int = cbstep/distrange*real_maxval

  public nchunk
  public distzero, distrange, distmin, distmax
  public corezero, corerange, coremin, coremax
  public lwpzero, lwprange, lwpmin, lwpmax
  public rwpzero, rwprange, rwpmin, rwpmax
  public thermzero, thermrange, thermmin, thermmax

  public i_corethres, i_lwpthres, i_thermthres, i_rwpthres

  ! tracking
  public n_minparentel

  ! splitting
  public n_growth_steps_min

  !main
  public rwpthres, lwpthres, corethres, thermthres

  public nmincells, nmincells_cloud

  public cbstep_as_int
  public real_maxval

end module constants
