module field_loader
  use modnetcdf, only: netcdfvar, check, inquire_ncvar, read_ncvar
  use tracking_common, only: tstart, nt, nx, ny
  use tracking_common, only: simulation_id
  use constants, only: nchunk

  use netcdf, only: nf90_open, nf90_nowrite, nf90_close

  implicit none

  ! the following arrays are for storing loaded fields and the mask in
  integer(kind=4), dimension(:,:,:), allocatable :: bool
  integer(kind=2), dimension(:,:,:,:), allocatable :: var
  integer, parameter :: var_ibase = 1, var_itop = 2  ! base and top values are always index 1 and 2 in `var`

  private
  integer :: finput, finput2
  real, dimension(:,:,:), allocatable :: readfield, readfield2

  public read_named_input
  public bool, var, var_ibase, var_itop

contains
  subroutine lookup_scaling_values(var_name, value_offset, value_scaling, min_value)
    use constants, only: distzero, distrange, distmin
    use constants, only: corezero, corerange, coremin
    use constants, only: lwpzero, lwprange, lwpmin
    use constants, only: rwpzero, rwprange, rwpmin

    character(len=*), intent(in) :: var_name

    real, intent(out) :: value_offset, value_scaling, min_value

    if (trim(var_name) == "core") then
      value_offset = corezero
      value_scaling = corerange
      min_value = coremin
    else if (trim(var_name) == "lwp") then
      value_offset = lwpzero
      value_scaling = lwprange
      min_value = lwpmin
    else if (trim(var_name) == "rwp") then
      value_offset = rwpzero
      value_scaling = rwprange
      min_value = rwpmin
    else if (trim(var_name) == "cldbase" .or. trim(var_name) == "cldtop") then
      value_offset = distzero
      value_scaling = distrange
      min_value = distmin
    else
      print *, "scaling values not defined for field `"//trim(var_name)//"`"
      call exit(3)
    endif
  end subroutine lookup_scaling_values

  !> Read field with name `nc_var%name` and rescale the field
  subroutine read_named_input(output, nc_var, nc_var2)
    use modnetcdf, only: fillvalue_i16

    integer(kind=2), dimension(:,:,:), intent(out) :: output
    type(netcdfvar), intent(in) :: nc_var
    type(netcdfvar), intent(in), optional :: nc_var2

    character(len=200) :: filename
    real :: value_offset, value_scaling, min_value
    integer :: k, kkmax, kk, j, i

    if (.not. allocated(readfield)) then
      allocate(readfield(nx, ny, nchunk))
    endif

    call lookup_scaling_values(nc_var%name, value_offset=value_offset, value_scaling=value_scaling,&
       min_value=min_value)

    filename = trim(simulation_id)//trim(nc_var%name)//'.nc'
    write(*,*) 'Reading ', trim(filename)
    call check ( nf90_open (trim(filename), NF90_NOWRITE, finput) )
    call inquire_ncvar(finput, nc_var)

    if (present(nc_var2)) then
      filename = trim(simulation_id)//trim(nc_var2%name)//'.nc'
      write(*,*) 'Reading second field ', trim(filename)
      call check ( nf90_open (trim(filename), NF90_NOWRITE, finput2) )
      call inquire_ncvar(finput2, nc_var2)
      if (.not. allocated(readfield2)) then
        allocate(readfield2(nx, ny, nchunk))
      endif
    endif

    do k = tstart,nt,nchunk
      kkmax = min(nt-k+1,nchunk)
      call read_ncvar(finput, nc_var, readfield,(/1,1,k/),(/nx,ny,kkmax/))
      if (present(nc_var2)) then
        call read_ncvar(finput2, nc_var2, readfield2,(/1,1,k/),(/nx,ny,kkmax/))
        readfield = readfield + readfield2
      endif

      do kk = 1,kkmax
        do j = 1, ny
          do i = 1, nx
            if (readfield(i,j,kk) >= min_value) then
              output(i,j,k+kk-1) = (readfield(i,j,kk) - value_offset ) / value_scaling
            else
              output(i,j,k+kk-1) = fillvalue_i16
            end if
          end do
        end do
      end do
    end do
    call check(nf90_close(finput))
  end subroutine read_named_input
end module field_loader
