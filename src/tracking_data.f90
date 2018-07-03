module tracking_data
  use tracking_common, only: celltype

  implicit none

  private

  type(celltype), pointer :: tracked_cores
  type(celltype), pointer :: tracked_clouds
  type(celltype), pointer :: tracked_rain
  type(celltype), pointer :: tracked_thermals

  integer :: icore = -1, ilwp = -1, irain = -1, ithermal = -1
  integer :: ncores, nclouds, nrains, nthermals
  integer :: nvar

  public tracked_cores
  public tracked_clouds
  public tracked_rain
  public tracked_thermals

  public icore, ithermal, ilwp, irain
  public ncores, nrains, nthermals, nclouds
  public nvar
end module tracking_data
