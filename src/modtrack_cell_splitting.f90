module modtrack_cell_splitting
  use tracking_common, only: ibase, itop, ivalue
  use tracking_common, only: celltype, cellptr

  implicit none

  public splitcell, find_split_regions

  private

  integer :: totcloudsystemnr = 0

  type :: cell_subregion_type
     integer(kind=2), allocatable, dimension(:,:) :: point_loc
     !> for each (i,j,t) entry in `points_loc` there is a index
     integer(kind=2), allocatable, dimension(:) :: parent_index
  end type cell_subregion_type

  contains


  !> Expand "active" regions while ensuring that the cloud-base height doesn't change rapidly
  !>
  !> Steps:
  !> 1. start from list of points (and parent indexes) which mark the region
  !> 2. use the obj_mask to grow
  !> 3. return the list points which are the union of the initial region and grown into region
  subroutine grow_split_regions(cell, var_base, &
        obj_mask, points, num_points, num_points_per_new_cell &
        )
    use constants, only: n_minparentel, n_growth_steps_min

    type(celltype), pointer, intent(in)   :: cell
    integer(kind=2), dimension(:,:,:), intent(in) :: var_base

    integer(kind=4), dimension(:,:,:), intent(inout) :: obj_mask
    integer, dimension(:,:), intent(inout) :: points
    integer, intent(inout) :: num_points
    integer, dimension(:), intent(inout) :: num_points_per_new_cell

    integer, allocatable, dimension(:)     :: maxiter
    integer :: iter = -1
    integer :: n = -1

    integer, allocatable, dimension(:,:)   :: boundary_points
    integer :: num_boundary_points = -1
    integer :: num_points_before_expansion = -1

    integer :: n_start = -1, n_end = -1

    ! keep track of how many points were part of cell before expanding it
    num_points_before_expansion = num_points

    ! there can be at maximum the same number of points in the boundary as the cell has points, this could probably be a lot less
    ! for optimisation...
    allocate(boundary_points(4,cell%n_points))
    ! start by including all points (including the interior) in the boundary points. On the first iteration of `findneighbours` we
    ! will only retain the boundary points from the first pass of growth
    boundary_points(:,:num_points) = points(:,:num_points)
    num_boundary_points = num_points

    ! calculate maximum number of growth steps for each parent
    allocate(maxiter(cell%nsplitters))
    maxiter = max( &
      num_points_per_new_cell(1:cell%nsplitters)/(20*n_minparentel), &
      n_growth_steps_min &
    )

    !print *, "maxiter=", maxiter
    !print *, "points=", points(4,num_points-10:num_points+1)

    do iter = 1, maxval(maxiter)
      !print *, "expansion iter", iter, num_boundary_points
      if (num_boundary_points == 0) then
         ! there are no more points matching the expansion criterion to grow into
         exit
      else
        ! we want to iterate over only the points that are in the boundary which are the last ones
        ! to have been added to the list
        n_start = num_points - num_boundary_points
        n_end = num_points

        ! reset the counter for how many points will be the next expanded boundary
        num_boundary_points = 0
        do n = n_start+1, n_end
          ! the 4th element of points maps to the index of the cell that is splitting a given point, use this to find the number of
          ! iterations to do for this splitter
          if (iter > maxiter(points(4,n))) then
             cycle
          endif
          call findneighbour(from_point=points(:,n), var_base=var_base, &
                             obj_mask=obj_mask, boundary_points=boundary_points, &
                             num_points_per_new_cell=num_points_per_new_cell, &
                             num_boundary_points=num_boundary_points)
        end do
      endif

      ! copy new boundary points onto list of points, i.e. the points that the original region has
      ! grown into
      points(:,num_points+1:num_points+num_boundary_points) = boundary_points(:,1:num_boundary_points)
      num_points = num_points + num_boundary_points
      !print *, "added", num_boundary_points
    end do

    !print *, "expanded from", num_points_before_expansion, "to", num_points, "points"
  end subroutine grow_split_regions


  !> For any points which were part of the oldcell but haven't been consumed by
  !> the regions-split-by-parents (or what these have expanded into) these points
  !> will be added to so-called "passive" cells
  !> The object mask is updated to indicate that these points have been added to
  !> cells
  subroutine create_cells_for_points_outside_expanded_regions(oldcell, obj_mask, &
         current_cell_points_loc, num_new_cells, &
         num_points_per_new_cell, new_points, num_new_points)

    use tracking_common, only: UNCLASSIFIED_IN_OBJECT

    type(celltype), intent(in), pointer     :: oldcell
    integer(kind=2), dimension(:,:), intent(in) :: current_cell_points_loc
    integer(kind=4), dimension(:,:,:), intent(inout) :: obj_mask
    integer, intent(inout) :: num_new_cells

    integer, intent(inout), dimension(:,:)  :: new_points
    integer, intent(inout) :: num_new_points
    integer, intent(inout), dimension(:)    :: num_points_per_new_cell

    integer :: n = -1
    integer(kind=2) :: i = -1, j = -1, t = -1

    integer, allocatable, dimension(:,:)   :: passive_points
    integer :: num_passive_points = 0

    integer, parameter :: MIN_NUM_POINTS_PASSIVE_CELL = 4

    ! only necessary to allocate for the number of points that haven't been added to "active" cells
    allocate(passive_points(4,oldcell%n_points-num_new_points))

    ! start off with as many cells as there were parent cells splitting the current cell we're
    ! working on
    num_new_cells = oldcell%nsplitters

    do n = 1, oldcell%n_points
      i = current_cell_points_loc(1,n)
      j = current_cell_points_loc(2,n)
      t = current_cell_points_loc(3,n)

      if (obj_mask(i,j,t) == UNCLASSIFIED_IN_OBJECT) then
        num_passive_points = 0
        num_new_cells = num_new_cells + 1

        call newpassive(i=i, j=j, t=t, &
          obj_mask=obj_mask, passive_points=passive_points, &
          num_passive_points=num_passive_points, num_new_cells=num_new_cells)

        if (num_passive_points > MIN_NUM_POINTS_PASSIVE_CELL) then
          new_points(:,num_new_points+1:num_new_points+num_passive_points) = passive_points(:,1:num_passive_points)
          num_points_per_new_cell(num_new_cells) = num_passive_points
          num_new_points = num_new_points + num_passive_points
        else
          num_new_cells = num_new_cells - 1
        endif
      end if
    end do

    !print *, "num passive cells added", num_new_cells - oldcell%nsplitters
  end subroutine


  !> Use a parent
  subroutine splitcell(cell, ncells, obj_mask, var_base, var_top, var_value, &
                       current_cell_points_loc, new_points, num_new_points, num_points_per_new_cell)
    use tracking_common, only: createcell
    use tracking_common, only: deletecell
    use tracking_common, only: print_cell_debug
    use tracking_common, only: UNCLASSIFIED_IN_OBJECT, PROCESSED_OBJECT
    use constants, only: n_minparentel

    type(celltype), pointer, intent(inout)   :: cell
    integer, intent(inout)           :: ncells
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
    integer, intent(inout), dimension(:,:)  :: new_points
    !> number of points which have been allocated into `new_points`
    integer, intent(inout) :: num_new_points
    integer, intent(inout), dimension(:)    :: num_points_per_new_cell

    type(celltype), pointer                  :: oldcell
    type(cellptr), allocatable, dimension(:) :: newcells
    integer :: n
    integer :: nn
    !> Counter for the number of cells we have after splitting
    integer :: num_new_cells

    integer, allocatable, dimension(:) :: cell_point_index_offset

    type(celltype), pointer :: parent_splitting_cell

    !print *, "before splitcell"
    !call print_cell_debug(cell)
    !if(cell%n_points> 10000000) write (*,*) 'Begin Splitcell'
    !write (*,*) 'Split cell with', cell%n_points, 'points in', cell%nsplitters, ' parts'

    oldcell => cell

    ! we start of with as many total cells as there are splitting cells, these will be the
    ! active ones
    num_new_cells = oldcell%nsplitters

    call grow_split_regions(cell=cell, var_base=var_base, &
                     obj_mask=obj_mask, points=new_points, num_points=num_new_points, &
                     num_points_per_new_cell=num_points_per_new_cell &
    )

    if (num_new_points < oldcell%n_points) then
      call create_cells_for_points_outside_expanded_regions(&
             oldcell=oldcell, obj_mask=obj_mask, &
             current_cell_points_loc=current_cell_points_loc, &
             num_new_cells=num_new_cells, &
             num_points_per_new_cell=num_points_per_new_cell, &
             new_points=new_points, num_new_points=num_new_points)
    end if

    ! allocatate storage for all new cells, active + passive
    allocate(newcells(num_new_cells))

    if (num_new_points > oldcell%n_points) then
       print *, "New cells from split cell have more points than original cell!"
       call exit(-1)
    endif

    ! allocate datastructures for storing new active cells
    do n = 1, oldcell%nsplitters
      !'Create and fill new active cell', n, '/', oldcell%nsplitters, "with", num_points_per_new_cell(n)

      ! this new active cell is split by one parent cell, which is the same parent cell as the old field was split by
      parent_splitting_cell => oldcell%splitters(n)%p

      call createcell(cell)
      newcells(n)%p => cell
      cell%nsplitters = 1
      allocate(cell%splitters(1))
      cell%splitters(1)%p => parent_splitting_cell
      cell%cloudtype = 4

      ! this parent splitting cell now only has one child cell
      parent_splitting_cell%nchildren = 1
      deallocate(parent_splitting_cell%children)
      allocate(parent_splitting_cell%children(1))
      parent_splitting_cell%children(1)%p => cell
    end do

    ! the passive cells will be the extra ones after we've added the number of "active" ones,
    ! add these now
    do n = oldcell%nsplitters + 1, num_new_cells
      call createcell(cell)
      newcells(n)%p => cell
      cell%cloudtype = 3
    end do

! 2018/07/09 (Leif): there is some old code here related to "siblings", unsure what this is. Seems incomplete
!
!       !Add the new cells to each other as siblings
!         do n = 1, num_new_cells
! !           newcells(n)%p%nsiblings = num_new_cells -1
!           allocate(newcells(n)%p%siblings(num_new_cells -1))
!           nnn = 0
!           do nn = 1, num_new_cells
!             if (n /= nn) then
!               nnn = nnn + 1
!               newcells(n)%p%siblings(nnn)%p => newcells(nn)%p
!             end if
!           end do
!         end do
!       end if

    !Allocate the arrays
    totcloudsystemnr = totcloudsystemnr + 1
    !print *, "creating output for a total of", num_new_cells, "new cells"
    do n = 1,num_new_cells
      !write (*,*) 'Create split cell nr ', n,'/',num_new_cells, 'with',num_points_per_new_cell(n),'elements'
      newcells(n)%p%n_points = num_points_per_new_cell(n)
      newcells(n)%p%cloudsystemnr = totcloudsystemnr
      allocate(newcells(n)%p%loc(3,num_points_per_new_cell(n)))
      allocate(newcells(n)%p%value(3,num_points_per_new_cell(n)))
    end do

    allocate(cell_point_index_offset(num_new_cells))
    cell_point_index_offset(:) = 0
    do n = 1, num_new_points
      !if (mod(n,1000000)==0) write (*,*) 'Fill new cells ', n,'/', num_new_points

      nn = new_points(4,n)
      cell_point_index_offset(nn) = cell_point_index_offset(nn) + 1
      newcells(nn)%p%loc(:,cell_point_index_offset(nn)) = new_points(1:3,n)
      newcells(nn)%p%value(ibase,cell_point_index_offset(nn)) = var_base(new_points(1,n),new_points(2,n),new_points(3,n))
      newcells(nn)%p%value(itop,cell_point_index_offset(nn)) = var_top(new_points(1,n),new_points(2,n),new_points(3,n))
      newcells(nn)%p%value(ivalue,cell_point_index_offset(nn)) = var_value(new_points(1,n),new_points(2,n),new_points(3,n))
    end do


    !Point at final new cell
    cell => newcells(num_new_cells)%p
    do n = 1,num_new_cells
      ncells = ncells + 1
      newcells(n)%p%id = ncells
    end do
    !Remove the original cell
    if(oldcell%n_points> 10000000) write (*,*) 'Split cell completed'
    call deletecell(oldcell)

    !print *, "after splitcell"
    !call print_cell_debug(cell)
  end subroutine splitcell


  !> Find the parent cell for each element in the current cell and update
  !> `splitters` and `children` in each. Also creates a number of each
  !> parent cell which assigned to the last row of `list`. `nr` appears to
  !> count the number of "child" elements that are split by each cell.
  !> 
  !> the cloud mask at given point gets set to the number of splitters
  subroutine find_split_regions(current_cell_points_loc, &
                                cell, parentarr, obj_mask, points, num_points, &
                                num_points_per_new_cell)

    use modarray, only: increase_array

    integer(kind=2), dimension(:,:), intent(in) :: current_cell_points_loc

    type(celltype),  pointer, intent(inout)                       :: cell
    type(cellptr),   allocatable, dimension(:,:,:), intent(inout) :: parentarr
    integer(kind=4), dimension(:,:,:), intent(inout)              :: obj_mask

    integer, dimension(:,:), intent(out)   :: points
    integer, dimension(:), intent(out)     :: num_points_per_new_cell
    integer, intent(out)                   :: num_points

    integer :: nn, i, j, t, n
    logical :: lnewparent

    do nn = 1, cell%n_points
      !write(*,*) 'Find splitter for element',nn,'/',cell%n_points
      if (mod(nn,1000000) == 0) write(*,*) 'Find splitter for element',nn,'/',cell%n_points
      i = current_cell_points_loc(1,nn)
      j = current_cell_points_loc(2,nn)
      t = current_cell_points_loc(3,nn)
      ! Is there a parent cell associated with this datapoint?
      if (associated(parentarr(i,j,t)%p)) then
        lnewparent = .true.
        ! Check if we've seen this cell before, i.e. has it be put into
        ! "cell.splitters" already?
        do n = 1, cell%nsplitters
          if (associated(cell%splitters(n)%p,parentarr(i,j,t)%p)) then
            lnewparent = .false.
            exit
          end if
        end do
        ! if it's new cell
        if (lnewparent) then
          ! add the parent to this cell's "splitters"
          cell%nsplitters = cell%nsplitters + 1
          call increase_array(cell%splitters, cell%nsplitters)
          cell%splitters(cell%nsplitters)%p => parentarr(i,j,t)%p

          ! and the current cell to "children" of the parent cell
          parentarr(i,j,t)%p%nchildren = parentarr(i,j,t)%p%nchildren + 1
          call increase_array(parentarr(i,j,t)%p%children, parentarr(i,j,t)%p%nchildren)
          parentarr(i,j,t)%p%children(parentarr(i,j,t)%p%nchildren)%p => cell

          n = cell%nsplitters
          ! the whole array is initiated to -1 to ensure we haven't added information for cells that
          ! don't exist, so set the counter to zero here
          num_points_per_new_cell(n) = 0
        end if

        num_points = num_points + 1
        obj_mask(i,j,t) = n
        points(1,num_points) = i
        points(2,num_points) = j
        points(3,num_points) = t
        points(4,num_points) = n
        num_points_per_new_cell(n) = num_points_per_new_cell(n) + 1
      end if
    end do
  end subroutine find_split_regions


  !> starting from position (in time and space) given as from (1,2,3) of the
  !> nnth element of `list` compare value at current position to neighbouring
  !> positions in time and space. The fourth component of `list` is used to
  !> assign into `obj_mask`
  subroutine findneighbour(from_point, var_base, obj_mask, boundary_points, &
                           num_points_per_new_cell, num_boundary_points)
    use constants, only: cbstep_as_int
    use tracking_common, only: tstart, nt, nx, ny
    use tracking_common, only: UNCLASSIFIED_IN_OBJECT

    integer,         intent(in), dimension(4)      :: from_point
    integer(kind=2), intent(in), dimension(:,:,:)    :: var_base

    integer(kind=4), intent(inout), dimension(:,:,:) :: obj_mask
    integer,         intent(inout), dimension(:,:)   :: boundary_points
    integer,         intent(inout), dimension(:)     :: num_points_per_new_cell
    integer,         intent(inout)                   :: num_boundary_points

    integer :: i, j, t
    integer :: ii ,jj, tt

    i = from_point(1)
    j = from_point(2)
    t = from_point(3)


    !Look west
    ii = i - 1
    if (ii.le.0) ii = nx
    if (obj_mask(ii,j,t)==UNCLASSIFIED_IN_OBJECT .and. var_base(ii,j,t) < var_base(i,j,t) + cbstep_as_int) then
      obj_mask(ii,j,t) = from_point(4)
      num_boundary_points = num_boundary_points + 1
      boundary_points(1:3,num_boundary_points) = (/ii,j,t/)
      boundary_points(4,num_boundary_points) = from_point(4)
      num_points_per_new_cell(from_point(4)) = num_points_per_new_cell(from_point(4)) + 1
    end if
    !Look east
    ii = i + 1
    if (ii.gt.nx) ii = 1
    if (obj_mask(ii,j,t)==UNCLASSIFIED_IN_OBJECT .and. var_base(ii,j,t) < var_base(i,j,t) + cbstep_as_int) then
      obj_mask(ii,j,t) = from_point(4)
      num_boundary_points = num_boundary_points + 1
      boundary_points(1:3,num_boundary_points) = (/ii,j,t/)
      boundary_points(4,num_boundary_points) = from_point(4)
      num_points_per_new_cell(from_point(4)) = num_points_per_new_cell(from_point(4)) + 1
    end if
    !Look north
    jj = j - 1
    if (jj.le.0) jj = ny
    if (obj_mask(i,jj,t)==UNCLASSIFIED_IN_OBJECT .and. var_base(i,jj,t) < var_base(i,j,t) + cbstep_as_int) then
      obj_mask(i,jj,t) = from_point(4)
      num_boundary_points = num_boundary_points + 1
      boundary_points(1:3,num_boundary_points) = (/i,jj,t/)
      boundary_points(4,num_boundary_points) = from_point(4)
      num_points_per_new_cell(from_point(4)) = num_points_per_new_cell(from_point(4)) + 1
    end if
    !Look south
    jj = j + 1
    if (jj.gt.ny) jj = 1
    if (obj_mask(i,jj,t)==UNCLASSIFIED_IN_OBJECT .and. var_base(i,jj,t) < var_base(i,j,t) + cbstep_as_int) then
      obj_mask(i,jj,t) = from_point(4)
      num_boundary_points = num_boundary_points + 1
      boundary_points(1:3,num_boundary_points) = (/i,jj,t/)
      boundary_points(4,num_boundary_points) = from_point(4)
      num_points_per_new_cell(from_point(4)) = num_points_per_new_cell(from_point(4)) + 1
    end if
    !Look forward
    tt = t+1
    if (tt <= nt) then
      if (obj_mask(i,j,tt)==UNCLASSIFIED_IN_OBJECT .and. var_base(i,j,tt) < var_base(i,j,t) + cbstep_as_int) then
        obj_mask(i,j,tt) = from_point(4)
        num_boundary_points = num_boundary_points + 1
        boundary_points(1:3,num_boundary_points) = (/i,j,tt/)
        boundary_points(4,num_boundary_points) = from_point(4)
        num_points_per_new_cell(from_point(4)) = num_points_per_new_cell(from_point(4)) + 1
      end if
    end if
    !Look backward
    tt = t - 1
    if (tt >=tstart) then
      if (obj_mask(i,j,tt)==UNCLASSIFIED_IN_OBJECT .and. var_base(i,j,tt) < var_base(i,j,t) + cbstep_as_int) then
        obj_mask(i,j,tt) = from_point(4)
        num_boundary_points = num_boundary_points + 1
        boundary_points(1:3,num_boundary_points) = (/i,j,tt/)
        boundary_points(4,num_boundary_points) = from_point(4)
        num_points_per_new_cell(from_point(4)) = num_points_per_new_cell(from_point(4)) + 1
      end if
    end if
  end subroutine findneighbour

  recursive subroutine newpassive(i, j, t, &
      obj_mask, passive_points, num_passive_points, num_new_cells)

    use tracking_common, only: tstart, nt, nx, ny
    use tracking_common, only: UNCLASSIFIED_IN_OBJECT, PROCESSED_OBJECT

    integer(kind=2), intent(in) :: i, j, t
    integer, intent(in) :: num_new_cells

    integer, intent(inout) :: num_passive_points
    integer, intent(inout), dimension(:,:)   :: passive_points
    integer(kind=4), intent(inout), dimension(:,:,:) :: obj_mask

    integer(kind=2) :: ii, jj, tt

    ! mark this point in the mast to indicate that it has been added to a cell (and so it won't be
    ! added to more than one cell)
    obj_mask(i,j,t) = PROCESSED_OBJECT

    ! add point to list of points for this passive cell
    num_passive_points = num_passive_points + 1
    passive_points(1:3,num_passive_points) = (/i,j,t/)
    passive_points(4,num_passive_points)   = num_new_cells

    !Look west
    ii = i - 1
    jj = j
    tt = t
    if (ii.le.0) ii = nx
    if (obj_mask(ii,jj,tt)==UNCLASSIFIED_IN_OBJECT) then
       call newpassive(i=ii, j=jj, t=tt, &
                       num_passive_points=num_passive_points, obj_mask=obj_mask, &
                       passive_points=passive_points, num_new_cells=num_new_cells)
    end if
    !Look east
    ii = i + 1
    jj = j
    tt = t
    if (ii.gt.nx) ii = 1
    if (obj_mask(ii,jj,tt)==UNCLASSIFIED_IN_OBJECT) then
       call newpassive(i=ii, j=jj, t=tt, &
                       num_passive_points=num_passive_points, obj_mask=obj_mask, &
                       passive_points=passive_points, num_new_cells=num_new_cells)
    end if

    !Look north
    ii = i
    jj = j - 1
    tt = t
    if (jj.le.0) jj = ny
    if (obj_mask(ii,jj,tt)==UNCLASSIFIED_IN_OBJECT) then
       call newpassive(i=ii, j=jj, t=tt, &
                       num_passive_points=num_passive_points, obj_mask=obj_mask, &
                       passive_points=passive_points, num_new_cells=num_new_cells)
    end if

    !Look south
    ii = i
    jj = j + 1
    tt = t
    if (jj.gt.ny) jj = 1
    if (obj_mask(ii,jj,tt)==UNCLASSIFIED_IN_OBJECT) then
       call newpassive(i=ii, j=jj, t=tt, &
                       num_passive_points=num_passive_points, obj_mask=obj_mask, &
                       passive_points=passive_points, num_new_cells=num_new_cells)
    end if

    !Look forward
    ii = i
    jj = j
    tt = t + 1
    if (tt <=nt) then
      if (obj_mask(ii,jj,tt)==UNCLASSIFIED_IN_OBJECT) then
       call newpassive(i=ii, j=jj, t=tt, &
                       num_passive_points=num_passive_points, obj_mask=obj_mask, &
                       passive_points=passive_points, num_new_cells=num_new_cells)
      end if
    end if

    !Look backward
    ii = i
    jj = j
    tt = t - 1
    if (tt>=tstart) then
      if (obj_mask(ii,jj,tt)==UNCLASSIFIED_IN_OBJECT) then
       call newpassive(i=ii, j=jj, t=tt, &
                       num_passive_points=num_passive_points, obj_mask=obj_mask, &
                       passive_points=passive_points, num_new_cells=num_new_cells)
      end if
    end if
  end subroutine newpassive
endmodule modtrack_cell_splitting
