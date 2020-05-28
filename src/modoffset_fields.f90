module offset_fields
  use tracking_common, only: dx, dy

  implicit none

  public

  !> offset velocities to correct for ambient wind or Galilean transform
  real :: u_vel_offset = 0, v_vel_offset = 0

  contains
    !> Offset the provided `input_field` by `time`, optionally getting back the
    !> original field from one already offset by passing in `add_offset=.true.`
    function offset_field(input_field, time, add_offset) result(output_field)
      real, dimension(:,:), intent(in) :: input_field
      real, intent(in) :: time
      logical, optional, intent(in) :: add_offset
      real, dimension(size(input_field,1),size(input_field,2)) :: output_field

      integer :: d, horz_index_shift

      horz_index_shift = 0
      d = 0
      output_field(:,:) = input_field(:,:)

      if (present(add_offset) .and. add_offset) then
        d = -1
      else
        d = 1
      endif

      if (u_vel_offset .ne. 0.0) then
        horz_index_shift = d*int(u_vel_offset*time/dx)
        output_field = cshift(output_field, shift=horz_index_shift, dim=1)
      endif

      if (v_vel_offset .ne. 0.0) then
        horz_index_shift = d*int(v_vel_offset*time/dy)
        output_field = cshift(output_field, shift=horz_index_shift, dim=2)
      endif
    end function offset_field

end module offset_fields
