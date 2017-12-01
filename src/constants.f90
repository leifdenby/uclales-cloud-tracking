module constants
  implicit none
  private

  real, parameter :: thermmin = -1.
  real, parameter :: thermmax = 10000.
  real, parameter :: thermthres = 300.

  real, parameter :: lwpmin   = -1.
  real, parameter :: lwpmax   = 10.
  real, parameter :: lwpthres = 0.01

  real, parameter :: coremin  = -5.
  real, parameter :: coremax  = 5.
  real, parameter :: corethres = 0.5

  real, parameter :: rwpmin   = -1.
  real, parameter :: rwpmax   = 10.
  real, parameter :: rwpthres = 0.01

  real, parameter :: distmin = -1.
  real, parameter :: distmax = 5000.

  real, parameter :: maxheight = 5000.

  integer, parameter :: nchunk = 100

  integer, parameter :: nmincells_cloud  = 1
  integer, parameter :: nmincells        = 4

  ! center and range values for parameters
  real, parameter :: thermzero  = 0.5*(thermmax + thermmin)
  real, parameter :: thermrange = (thermmax - thermmin)/real(huge(1_2))

  real, parameter :: lwpzero  = 0.5*(lwpmax + lwpmin)
  real, parameter :: lwprange = (lwpmax - lwpmin)/real(huge(1_2))

  real, parameter :: corezero  = 0.5*(coremax + coremin)
  real, parameter :: corerange = (coremax - coremin)/real(huge(1_2))

  real, parameter :: rwpzero  = 0.5*(rwpmax + rwpmin)
  real, parameter :: rwprange = (rwpmax - rwpmin)/real(huge(1_2))

  real, parameter :: distzero  = 0.5*(distmax + distmin)
  real, parameter :: distrange = (distmax - distmin)/real(huge(1_2))

  integer(kind=2), parameter :: i_lwpthres   = (lwpthres - lwpzero)/lwprange
  integer(kind=2), parameter :: i_corethres  = (corethres - corezero)/corerange
  integer(kind=2), parameter :: i_thermthres = (thermthres - thermzero)/thermrange
  integer(kind=2), parameter :: i_rwpthres   = (rwpthres - rwpzero)/rwprange

  public nchunk
  public distzero, distrange, distmin
  public corezero, corerange, coremin
  public lwpzero, lwprange, lwpmin
  public rwpzero, rwprange, rwpmin

end module constants
