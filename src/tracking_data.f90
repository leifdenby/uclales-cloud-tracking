module tracking_data
  use tracking_common, only: celltype

  implicit none

  private

  type(celltype), pointer :: core, cloud, rain, thermal

  integer :: icore = -1, ilwp = -1, irain = -1, ithermal = -1
  integer :: ncores, nclouds, nrains, nthermals
  integer :: nvar

  public core, cloud, rain, thermal

  public icore, ithermal, ilwp, irain
  public ncores, nrains, nthermals, nclouds
  public nvar
end module tracking_data
