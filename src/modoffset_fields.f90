module offset_fields
  use tracking_common, only: dx, dy

  implicit none

  public

  !> offset velocities to correct for ambient wind or Galilean transform
  real :: u_vel_offset = 0, v_vel_offset = 0

  interface offset_field
    module procedure offset_field_i
    module procedure offset_field_r
  end interface offset_field


  contains
    !> Offset the provided real `input_field` by `time`, optionally getting back the
    !> original field from one already offset by passing in `remove_offset=.true.`
    function offset_field_r(input_field, time, remove_offset) result(output_field)
      real, dimension(:,:), intent(in) :: input_field
      real, intent(in) :: time
      logical, optional, intent(in) :: remove_offset
      real, dimension(size(input_field,1),size(input_field,2)) :: output_field

      integer :: horz_index_shift

      horz_index_shift = 0
      output_field(:,:) = input_field(:,:)

      if (u_vel_offset .ne. 0.0) then
        horz_index_shift = calc_index_offset(time=time, offset_dim=1, remove_offset=remove_offset)
        output_field = cshift(output_field, shift=horz_index_shift, dim=1)
      endif

      if (v_vel_offset .ne. 0.0) then
        horz_index_shift = calc_index_offset(time=time, offset_dim=2, remove_offset=remove_offset)
        output_field = cshift(output_field, shift=horz_index_shift, dim=2)
      endif
    end function offset_field_r


    !> Offset the provided integer `input_field` by `time`, optionally getting back the
    !> original field from one already offset by passing in `remove_offset=.true.`
    function offset_field_i(input_field, time, remove_offset) result(output_field)
      integer, dimension(:,:), intent(in) :: input_field
      real, intent(in) :: time
      logical, optional, intent(in) :: remove_offset
      integer, dimension(size(input_field,1),size(input_field,2)) :: output_field

      integer :: horz_index_shift

      horz_index_shift = 0
      output_field(:,:) = input_field(:,:)

      if (u_vel_offset .ne. 0.0) then
        horz_index_shift = calc_index_offset(time=time, offset_dim=1, remove_offset=remove_offset)
        output_field = cshift(output_field, shift=horz_index_shift, dim=1)
      endif

      if (v_vel_offset .ne. 0.0) then
        horz_index_shift = calc_index_offset(time=time, offset_dim=2, remove_offset=remove_offset)
        output_field = cshift(output_field, shift=horz_index_shift, dim=2)
      endif
    end function offset_field_i


    function calc_index_offset(time, offset_dim, remove_offset) result(horz_index_shift)
      integer, intent(in) :: offset_dim
      real, intent(in) :: time
      logical, optional :: remove_offset

      integer :: d, horz_index_shift

      horz_index_shift = 0
      d = 0

      if (present(remove_offset) .and. remove_offset) then
        d = -1
      else
        d = 1
      endif

      if (offset_dim == 1) then
        horz_index_shift = d*int(u_vel_offset*time/dx)
      else if (offset_dim == 2) then
        horz_index_shift = d*int(v_vel_offset*time/dy)
      else
        print *, "Requested horz index shift for invalid dimension=", offset_dim
        call exit(2)
      endif
    end function calc_index_offset

end module offset_fields
