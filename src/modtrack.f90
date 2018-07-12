module modtrack
  use tracking_common, only: celltype, cellptr

  use tracking_common, only: dt, dx, dy
  use tracking_common, only: ibase, itop, ivalue
  use tracking_common, only: nt, nx, ny
  use tracking_common, only: tstart
  use tracking_common, only: createcell, deletecell, firstcell, nextcell
  use tracking_common, only: INSIDE_OBJECTS

  use constants, only: z_offset => parent_array_allowed_height_offset

  use modtrack_cell_splitting, only: splitcell

  use modarray, only: increase_array

  implicit none

  private

  ! main
  public dotracking, fillparentarr, findparents

  ! modstatistics
  public nextcell, firstcell

  contains

  subroutine findparents(cell, parentarr, base, top)
    use constants, only: distrange, real_maxval

    type(celltype), pointer, intent(inout)         :: cell
    type(cellptr), allocatable,dimension(:,:,:), intent(inout) :: parentarr
    integer(kind=2), dimension(:,:,:),allocatable, intent(in)             :: base, top
    integer :: n, nn, i, j, t, iret, count
    logical :: lnewparent

    count = 0
    write (*,*) '.. entering findparents'
    iret = firstcell(cell)
    do
      if (iret == -1) exit
      count = count +1
      do nn = 1, cell%n_points
        i = cell%loc(1,nn)
        j = cell%loc(2,nn)
        t = cell%loc(3,nn)
        if (associated(parentarr(i,j,t)%p)) then
          if (& ! parent cell's base must be below top of cell, avoid floating thermals above cloud
              base(i,j,t) <= cell%value(itop,nn)  .and. &
              ! top of cell must be above base of cloud, this is the requirement of overlap
              top(i,j,t) + z_offset*real_maxval/distrange >= cell%value(ibase,nn) &
                 ) then
            lnewparent = .true.
            do n = 1, cell%nparents
              if (cell%parents(n)%p%id == parentarr(i,j,t)%p%id) then
                lnewparent = .false.
                exit
              end if
            end do
            if (lnewparent) then
              cell%nparents = cell%nparents + 1
              call increase_array(cell%parents, cell%nparents)
              cell%parents(cell%nparents)%p => parentarr(i,j,t)%p

              parentarr(i,j,t)%p%nchildren = parentarr(i,j,t)%p%nchildren + 1
              call increase_array(parentarr(i,j,t)%p%children, parentarr(i,j,t)%p%nchildren)
              parentarr(i,j,t)%p%children(parentarr(i,j,t)%p%nchildren)%p => cell
            end if
          end if
        end if
      end do
      iret = nextcell(cell)
    end do
    write (*,*) '.. leaving findparents'
  end subroutine findparents

  !> Iterate over all data-points in space and time and construct `celltype`
  !> instances (stored in a linked list through the `next`/`previous` attributes) 
  subroutine dotracking(cell, ncells, nmincells, obj_mask, var_base, var_top, var_value, parentarr)
    type(celltype), pointer, intent(inout)                   :: cell
    integer, intent(out)                             :: ncells
    integer, intent(in)                                      :: nmincells
    integer(kind=2), dimension(:,:,:), intent(in) :: var_base
    integer(kind=2), dimension(:,:,:), intent(in) :: var_top
    integer(kind=2), dimension(:,:,:), intent(in) :: var_value

    type(cellptr), allocatable, dimension(:,:,:), intent(inout), optional :: parentarr
    integer(kind=4), dimension(:,:,:), intent(inout) :: obj_mask

    integer(kind=2), dimension(:,:), allocatable :: current_cell_points_loc

    integer :: i, j, t
    write (*,*) '.. entering tracking'

    allocate(current_cell_points_loc(3,ceiling(min(0.3*(huge(1)-2),0.5*real(nx)*real(ny)*real(nt-tstart)))))
    print *, "Allocating array for storing cell locations"

    nullify(cell)
    do t = tstart, nt
      if(mod(t,10)==0) write (*,'(A,I10,A,I10)') "Time =", t,"  Nr of Cells", ncells
      do  j = 1, ny
        do i = 1, nx
          if (obj_mask(i,j,t)==INSIDE_OBJECTS) then
             !write (*,*) 'Create a new cell'
            call createcell(cell)
            call identify_new_cell(i, j, t, cell, obj_mask, var_base, var_top, current_cell_points_loc=current_cell_points_loc)
            if (cell%n_points >= nmincells) then
              if(cell%n_points> 10000000) write (*,*) '..finalizing cell'
              call finalizecell(cell, ncells, parentarr, obj_mask, var_base, var_top, &
                                var_value, current_cell_points_loc=current_cell_points_loc)
            else
              call deletecell(cell)
            end if
          end if
        end do
      end do
    end do
    write (*,*) '.. leaving tracking'
    deallocate(current_cell_points_loc)
  end subroutine dotracking


  !> Copy the information from the temporary array `current_cell_points_loc` into the `cell`
  !> instance's `loc` (location) and `value` attributes.
  !!
  !! - `current_cell_points_loc` contains all the positions in space of time of datapoints which
  !! are part of the current cell
  !! - copy into `cell.loc` the location in space of time of each element
  !! - set in `cell.value` the scalar values of each element (which a defined in
  !! both space and time) wrt cloud-top, cloud-base and "ivalue" ??
  !! - obj_mask is set to `-1` for all datapoints for the elements in the current
  !! cell
  !!
  !! @TODO what does `ivalue` mean here?
  subroutine finalizecell(cell, ncells, parentarr, obj_mask, var_base, var_top, var_value, current_cell_points_loc)
    use tracking_common, only: PROCESSED_OBJECT
    use modtrack_cell_splitting, only: find_split_regions

    type(celltype), pointer, intent(inout)                   :: cell
    integer, intent(out)                             :: ncells
    type(cellptr), allocatable, dimension(:,:,:), intent(inout), optional :: parentarr
    integer :: n
    integer(kind=4), dimension(:,:,:), intent(inout) :: obj_mask
    integer(kind=2), dimension(:,:,:), intent(in) :: var_base
    integer(kind=2), dimension(:,:,:), intent(in) :: var_top
    integer(kind=2), dimension(:,:,:), intent(in) :: var_value
    integer(kind=2), dimension(:,:), intent(in) :: current_cell_points_loc

    !> for storing the points of the new split regions (including points grown into by expand_regions
    !> below)
    !> list has shape (4, cell%n_points) and stores (i,j,t,n) where i,j and the
    !> position indecies, t the time index and n the number of splitters for this
    !> element of the cell
    integer, allocatable, dimension(:,:)  :: new_points
    !> number of points which have been allocated into `new_points`
    integer :: num_new_points
    integer, allocatable, dimension(:)    :: num_points_per_new_cell

    logical :: cell_was_split

    integer :: i, j, tn

    ! allocate storage for new points, for storing how many points each parent cell will have and
    ! counter to indicate how much of storage has been used
    allocate(new_points(4,cell%n_points))
    allocate(num_points_per_new_cell(cell%n_points/2))
    num_new_points = 0
    new_points(:,:) = -1
    num_points_per_new_cell(:) = -1

    cell_was_split = .false.


    if (present(parentarr)) then
      ! find out how the cells in the parent dataset overlap with the current cell to determine
      ! whether we will want to split this cell, the subroutine below also returns data about the
      ! split regions so these can be use when creating the split-off cells
      call find_split_regions(current_cell_points_loc=current_cell_points_loc, &
                              cell=cell, parentarr=parentarr, obj_mask=obj_mask, &
                              points=new_points, num_points=num_new_points, &
                              num_points_per_new_cell=num_points_per_new_cell  &
                             )
      if (cell%nsplitters >= 2) then
        call splitcell(cell, ncells, obj_mask, var_base, var_top, var_value, &
                       current_cell_points_loc=current_cell_points_loc, new_points=new_points, &
                       num_new_points=num_new_points, num_points_per_new_cell=num_points_per_new_cell)

        cell_was_split = .true.
      elseif (cell%nsplitters == 1) then
       cell%cloudtype = 1
      elseif (cell%nsplitters == 0) then
       cell%cloudtype = 2
      else
       print *, "Error: not implemented"
       call exit(-2)
      endif
    end if

    if (.not. cell_was_split) then
      allocate(cell%loc(3,cell%n_points))
      allocate(cell%value(3,cell%n_points))
      cell%loc(:,1:cell%n_points) = current_cell_points_loc(:,1:cell%n_points)

      do n = 1, cell%n_points
        i = current_cell_points_loc(1,n)
        j = current_cell_points_loc(2,n)
        tn = current_cell_points_loc(3,n)
        cell%value(ibase, n) = var_base(i, j, tn)
        cell%value(itop, n)  = var_top(i, j, tn)
        cell%value(ivalue, n) = var_value(i, j, tn)
        obj_mask(i, j, tn) = PROCESSED_OBJECT
      end do

      ncells = ncells + 1
      cell%id = ncells
    endif
  end subroutine finalizecell

  !> Given the current location (i,j,t) in space and time look at neighbouring
  !> points in space and time and if they satisfy the constraints from being part of
  !> the same cell:
  !!
  !! set obj_mask=-2 so that this data-point is not considered twice for multiple
  !! cells, store the position in space and time into the `current_cell_points_loc` array and
  !! increment the `n_points` counter on the provided cell
  !!
  !! TODO: Looking west/east/north/south and truncating the indexing near the
  !! edge is unnecessarily costly, the current element will be checked twice
  recursive subroutine identify_new_cell(i, j, t, cell, obj_mask, var_base, var_top, current_cell_points_loc)
    use tracking_common, only: UNCLASSIFIED_IN_OBJECT, INSIDE_OBJECTS

    integer, intent(in)                                      :: i, j, t
    type(celltype),pointer, intent(inout)                    :: cell
    integer :: ii, jj, tt
    integer(kind=4), dimension(:,:,:), intent(inout) :: obj_mask
    integer(kind=2), dimension(:,:), intent(inout) :: current_cell_points_loc

    integer(kind=2), dimension(:,:,:), intent(in) :: var_base
    integer(kind=2), dimension(:,:,:), intent(in) :: var_top

    cell%n_points = cell%n_points + 1
    if (mod(cell%n_points,10000000) == 0) write(*,*) 'Cell element ', cell%n_points
    current_cell_points_loc(1,cell%n_points) = i
    current_cell_points_loc(2,cell%n_points) = j
    current_cell_points_loc(3,cell%n_points) = t

    obj_mask(i, j, t) = UNCLASSIFIED_IN_OBJECT

    !Look west
    ii = i - 1
    jj = j
    tt = t
    if (ii.le.0) ii = nx
    if (obj_mask(ii,jj,tt)==INSIDE_OBJECTS) then
      if (var_base(i,j,t) <= var_top(ii,jj,tt) .and. var_top(i,j,t) >= var_base(ii,jj,tt)) then
        call identify_new_cell(ii,jj,tt,cell, obj_mask, var_base, var_top, current_cell_points_loc=current_cell_points_loc)
      end if
    end if
    !Look east
    ii = i + 1
    jj = j
    tt = t
    if (ii.gt.nx) ii = 1
    if (obj_mask(ii,jj,tt)==INSIDE_OBJECTS) then
      if (var_base(i,j,t) <= var_top(ii,jj,tt) .and. var_top(i,j,t) >= var_base(ii,jj,tt)) then
        call identify_new_cell(ii,jj,tt,cell, obj_mask, var_base, var_top, current_cell_points_loc=current_cell_points_loc)
      end if
    end if

    !Look north
    ii = i
    jj = j - 1
    tt = t
    if (jj.le.0) jj = ny
    if (obj_mask(ii,jj,tt)==INSIDE_OBJECTS) then
      if (var_base(i,j,t) <= var_top(ii,jj,tt) .and. var_top(i,j,t) >= var_base(ii,jj,tt)) then
        call identify_new_cell(ii,jj,tt,cell, obj_mask, var_base, var_top, current_cell_points_loc=current_cell_points_loc)
      end if
    end if

    !Look south
    ii = i
    jj = j + 1
    tt = t
    if (jj.gt.ny) jj = 1
    if (obj_mask(ii,jj,tt)==INSIDE_OBJECTS) then
      if (var_base(i,j,t) <= var_top(ii,jj,tt) .and. var_top(i,j,t) >= var_base(ii,jj,tt)) then
        call identify_new_cell(ii,jj,tt,cell, obj_mask, var_base, var_top, current_cell_points_loc=current_cell_points_loc)
      end if
    end if

    !Look forward
    ii = i
    jj = j
    tt =t + 1
    if (tt <= nt) then
      if (obj_mask(ii,jj,tt)==INSIDE_OBJECTS) then
        if (var_base(i,j,t) <= var_top(ii,jj,tt) .and. var_top(i,j,t) >= var_base(ii,jj,tt)) then
          call identify_new_cell(ii,jj,tt,cell, obj_mask, var_base, var_top, current_cell_points_loc=current_cell_points_loc)
        end if
      end if
    end if
    !Look backward
    ii = i
    jj = j
    tt = t - 1
    if (tt >= tstart) then
      if (obj_mask(ii,jj,tt)==INSIDE_OBJECTS) then
        if (var_base(i,j,t) <= var_top(ii,jj,tt) .and. var_top(i,j,t) >= var_base(ii,jj,tt)) then
          call identify_new_cell(ii,jj,tt,cell, obj_mask, var_base, var_top, current_cell_points_loc=current_cell_points_loc)
        end if
      end if
    end if
  end subroutine identify_new_cell

  !> Create a 3D array mapping from all datapoints in space and time to the cell
  !>  associated with each datapoint, but only for cells where the cell's base
  !> height is below what is passed in as `minbase`. This will for example
  !> provide a mapping from (x,y,t) to the "core" cells where the core cell
  !> reaches below a given height
  !>
  !! Also optionally store the "cloud-base" and "cloud-top" value into `base` and
  !! `top` arrays
  !!
  !! @TODO Does this subroutine actually read from `parentarr`?
  subroutine fillparentarr(cell, minbase, parentarr, base, top)
    use modnetcdf, only : fillvalue_i16
    use constants, only : n_minparentel

    type(celltype), pointer, intent(inout)         :: cell
    integer(kind=2), dimension(:), intent(in)      :: minbase
    type(cellptr), dimension(:,:,:), intent(inout)   :: parentarr
    integer(kind=2), dimension(:,:,:), intent(out), optional  :: base, top

    integer :: c = 0

    integer :: n, i, j, t, iret
    write (*,*) '.. entering fillparentarr'
    if (present(base)) then
      base = fillvalue_i16
      top  = fillvalue_i16
    end if
    do t = tstart, nt
      do j = 1, ny
        do i = 1, nx
          nullify (parentarr(i,j,t)%p)
        end do
      end do
    end do
    iret = firstcell(cell)
    do
      if (iret == -1) exit
      if (cell%n_points > n_minparentel) then
        if (minval(cell%value(ibase,1:cell%n_points))< minbase(cell%loc(3,1))) then
          do n = 1, cell%n_points
            parentarr(cell%loc(1,n),cell%loc(2,n),cell%loc(3,n))%p => cell
            if (present(base)) then
              base(cell%loc(1,n),cell%loc(2,n),cell%loc(3,n))        =  cell%value(ibase,n)
              top(cell%loc(1,n),cell%loc(2,n),cell%loc(3,n))        =  cell%value(itop,n)
            end if
            c = c+1
          end do
        end if
      end if
      iret =  nextcell(cell)
    end do
    write (*,*) '.. leaving fillparentarr', c
  end subroutine fillparentarr


  !> Check if successive timesteps appear to strange number of `obj_mask==0`
  !> data-points
  !!
  !! Strange defined as either:
  !! - number of cells with obj_mask==0 was above 100 in previous step and is 0 in current 
  !! - succesive timesteps has less than 25% of obj_mask==0 than previous timestep
  !!
  !! There's something strange whereby obj_mask is actually changed if either of the
  !! above is true, I think this is to always relate to the last timestep that
  !! was deemed to be good for comparison, maybe...
  subroutine checkframes(obj_mask)
    use tracking_common, only: INSIDE_OBJECTS

    integer(kind=4), dimension(:,:,:), intent(inout) :: obj_mask

    integer :: t, count1, count2
    real    :: hlp

      !..check for empty frames
    write(0,*) 'Check frames'
    count2 = 0
    do t = 2, nt
      count1 = count(obj_mask(:,:,t)==INSIDE_OBJECTS)
      hlp = float(count1)/float(count2)
      if (count1.eq.0.and.count2.gt.100) then
          write (0,*) 'WARNING: empty frame at t = ',t
          obj_mask(:,:,t)= obj_mask(:,:,t-1)
          count1 = count2
      elseif (hlp.lt.0.25) then
          write (0,*) 'WARNING: strange frame at t = ',t
          obj_mask(:,:,t)= obj_mask(:,:,t-1)
          count1 = count2
      end if
      count2 = count1
    end do

  end subroutine checkframes
end module modtrack
