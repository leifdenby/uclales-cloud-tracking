module modtrack_cell_splitting
  use tracking_common, only: ibase, itop, ivalue
  use tracking_common, only: celltype, cellptr

  implicit none

  public splitcell

  private

  integer :: totcloudsystemnr = 0

  type :: cell_subregion_type
     integer(kind=2), allocatable, dimension(:,:) :: point_loc
     !> for each (i,j,t) entry in `points_loc` there is a index
     integer(kind=2), allocatable, dimension(:) :: parent_index
  end type cell_subregion_type

  contains

  !> Find the parent cell for each element in the current cell and update
  !> `splitters` and `children` in each. Also creates a number of each
  !> parent cell which assigned to the last row of `list`. `nr` appears to
  !> count the number of "child" elements that are split by each cell.
  !> 
  !> the cloud mask at given point gets set to the number of splitters
  subroutine identify_splitting_overlap(cell, parentarr, &
                                 points, &
                                 num_points_per_new_cell, n_split_points, &
                                 current_cell_points_loc)

    use modarray, only: increase_array

    type(celltype),  pointer, intent(inout)                       :: cell
    type(cellptr),   allocatable, dimension(:,:,:), intent(inout) :: parentarr
    integer(kind=2), dimension(:,:), intent(in) :: current_cell_points_loc

    integer, dimension(:,:), intent(out)   :: points
    integer, dimension(:), intent(out)     :: num_points_per_new_cell
    integer, intent(inout)                   :: n_split_points

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
        end if
        n_split_points = n_split_points + 1
        points(1,n_split_points) = i
        points(2,n_split_points) = j
        points(3,n_split_points) = t
        points(4,n_split_points) = n
        num_points_per_new_cell(n) = num_points_per_new_cell(n) + 1
      end if
    end do
  end subroutine identify_splitting_overlap


  subroutine grow_split_regions(cell, var_base, &
        obj_mask, points, num_points, num_points_per_new_cell &
        )
    ! 1. start from list of points (and parent indexes) which mark the region
    ! 2. use the obj_mask to grow
    ! 3. return the list points which are the union of the initial region and grown into region
    use constants, only: n_minparentel, n_growth_steps_max

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
      num_points_per_new_cell(1:cell%nsplitters)/(10*n_minparentel), &
      n_growth_steps_max &
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
             call exit(-1)
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
  !subroutine create_cells_for_points_outside_expanded_regions(oldcell, current_cell_points_loc, &
         !obj_mask, num_points_per_new_cell)

    !use tracking_common, only: UNCLASSIFIED_IN_OBJECT

    !integer(kind=2), dimension(:,:), intent(in) :: current_cell_points_loc
    !type(celltype), pointer, intent(in)   :: oldcell
    !integer(kind=4), dimension(:,:,:), intent(inout) :: obj_mask
    !integer, dimension(:), intent(inout)    :: num_points_per_new_cell

    !integer :: n = -1
    !integer :: i, j, t

    !integer :: totnewcells = -1
    !integer :: nractive = -1

    !loop: do n = 1, oldcell%n_points
      !i = current_cell_points_loc(1,n)
      !j = current_cell_points_loc(2,n)
      !t = current_cell_points_loc(3,n)
      !if (obj_mask(i,j,t) == UNCLASSIFIED_IN_OBJECT) &
      !then !Found new outflow region

!!         npassive = npassive + 1
        !totnewcells = totnewcells + 1
        !num_points_per_new_cell(totnewcells) = 0
        !nractive = -1
        !call newpassive(i=i, j=j, t=t, nr=num_points_per_new_cell, totnewcells=totnewcells, &
                        !nendlist=nendlist, obj_mask=obj_mask, endlist=endlist, nractive=nractive)
        !if (nractive > 0) then
          !endlist(4,1+nendlist-nr(totnewcells):nendlist) = nractive
          !num_points_per_new_cell(nractive) = num_points_per_new_cell(nractive) + num_points_per_new_cell(totnewcells)
          !nr(totnewcells) = 0
          !totnewcells = totnewcells -1
        !end if
        !if (mod(n,print_interval)==0) write (*,*) 'Passive cell',n,'/',oldcell%n_points
      !end if
    !end do loop
  !end subroutine

  !> Use a parent
  subroutine splitcell(cell, ncells, parentarr, obj_mask, var_base, var_top, var_value, current_cell_points_loc)
    use tracking_common, only: nrel_max
    use tracking_common, only: createcell
    use tracking_common, only: deletecell
    use tracking_common, only: print_cell_debug
    use tracking_common, only: UNCLASSIFIED_IN_OBJECT, PROCESSED_OBJECT
    use constants, only: n_minparentel

    type(celltype), pointer, intent(inout)   :: cell
    integer, intent(inout)           :: ncells
    type(cellptr), allocatable, dimension(:,:,:), intent(inout) :: parentarr
    integer(kind=4), dimension(:,:,:), intent(inout) :: obj_mask

    integer(kind=2), dimension(:,:,:), intent(in) :: var_base
    integer(kind=2), dimension(:,:,:), intent(in) :: var_top
    integer(kind=2), dimension(:,:,:), intent(in) :: var_value
    integer(kind=2), dimension(:,:), intent(in) :: current_cell_points_loc

    type(celltype), pointer                  :: oldcell
    type(cellptr), allocatable, dimension(:) :: newcells
    integer :: n, tmin, tmax
    integer :: nn, num_points, npassive, totnewcells

    !> for storing the points of the new split regions (including points grown into by expand_regions
    !> below)
    !> list has shape (4, cell%n_points) and stores (i,j,t,n) where i,j and the
    !> position indecies, t the time index and n the number of splitters for this
    !> element of the cell
    integer, allocatable, dimension(:,:)  :: new_points
    !> number of points which have been allocated into `new_points`
    integer :: num_new_points

    integer, allocatable, dimension(:)    :: num_points_per_new_cell

    integer :: print_interval = 1

    type(celltype), pointer :: parent_splitting_cell

    !print *, "before splitcell"
    !call print_cell_debug(cell)

    oldcell => cell
    if(cell%n_points> 10000000) write (*,*) 'Begin Splitcell'

    ! allocate storage for new points, for storing how many points each parent cell will have and
    ! counter to indicate how much of storage has been used
    allocate(new_points(4,cell%n_points))
    allocate(num_points_per_new_cell(cell%n_points/2))
    num_new_points = 0
    new_points(:,:) = 0
    num_points_per_new_cell(:) = 0

    num_points = 0
    npassive = 0
    totnewcells = 0
    !nr = 0

    tmin = minval(current_cell_points_loc(3,1:cell%n_points))
    tmax = maxval(current_cell_points_loc(3,1:cell%n_points))

    call find_split_regions(current_cell_points_loc=current_cell_points_loc, &
                            cell=cell, parentarr=parentarr, obj_mask=obj_mask, &
                            points=new_points, num_points=num_new_points, &
                            num_points_per_new_cell=num_points_per_new_cell  &
                           )

    if (oldcell%nsplitters >= 2) then
          write (*,*) 'Split cell with', cell%n_points, 'points in', cell%nsplitters, ' parts'

          call grow_split_regions(cell=cell, var_base=var_base, &
                           obj_mask=obj_mask, points=new_points, num_points=num_new_points, &
                           num_points_per_new_cell=num_points_per_new_cell &
          )

          !Add cells for the number of parentcells
      !       npassive = 0
          totnewcells = oldcell%nsplitters
          !if (num_new_points < oldcell%n_points) then
          !end if

          ! allocatate storage for all new cells, active + passive
          allocate(newcells(totnewcells))

          ! allocate datastructures for storing new active cells
          do n = 1, oldcell%nsplitters
            if(mod(n,print_interval)==0) write (*,*) 'Create and fill new active cell', n, '/', oldcell%nsplitters

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

          ! XXX: disable allocation of passive cells for now, REALLY needs rewriting
          !do n = oldcell%nsplitters + 1, totnewcells
            !if(mod(n,print_interval)==0) write (*,*) 'Create and fill new passive cell', n, '/', totnewcells
            !call createcell(cell)
            !newcells(n)%p => cell
            !cell%cloudtype = 3
          !end do
      !       if (lsiblings) then
      !             nrel_max = max(nrel_max, totnewcells)
      !
      !       !Add the new cells to each other as siblings
      !         do n = 1, totnewcells
      ! !           newcells(n)%p%nsiblings = totnewcells -1
      !           allocate(newcells(n)%p%siblings(totnewcells -1))
      !           nnn = 0
      !           do nn = 1, totnewcells
      !             if (n /= nn) then
      !               nnn = nnn + 1
      !               newcells(n)%p%siblings(nnn)%p => newcells(nn)%p
      !             end if
      !           end do
      !         end do
      !       end if
          !Allocate the arrays
          totcloudsystemnr = totcloudsystemnr + 1
          print *, "creating output for a total of", totnewcells, "new cells"
          do n = 1,totnewcells
            write (*,*) 'Create split cell nr ', n,'/',totnewcells, 'with',num_points_per_new_cell(n),'elements'
      !         currid = oldcell%splitters(n)%p%id
            newcells(n)%p%n_points = num_points_per_new_cell(n)
            newcells(n)%p%cloudsystemnr = totcloudsystemnr
            allocate(newcells(n)%p%loc(3,num_points_per_new_cell(n)))
            allocate(newcells(n)%p%value(3,num_points_per_new_cell(n)))
          end do
          num_points_per_new_cell = 0
          do n = 1, num_new_points
            if (mod(n,1000000)==0) write (*,*) 'Fill new cells ', n,'/', num_new_points
            nn = new_points(4,n)
            num_points_per_new_cell(nn) = num_points_per_new_cell(nn) + 1
            newcells(nn)%p%loc(:,num_points_per_new_cell(nn)) = new_points(1:3,n)
            newcells(nn)%p%value(ibase,num_points_per_new_cell(nn)) = var_base(new_points(1,n),new_points(2,n),new_points(3,n))
            newcells(nn)%p%value(itop,num_points_per_new_cell(nn)) = var_top(new_points(1,n),new_points(2,n),new_points(3,n))
            newcells(nn)%p%value(ivalue,num_points_per_new_cell(nn)) = var_value(new_points(1,n),new_points(2,n),new_points(3,n))
          end do


          !Point at final new cell
          cell => newcells(totnewcells)%p
          do n = 1,totnewcells - 1
            ncells = ncells + 1
            newcells(n)%p%id = ncells
          end do
          !Remove the original cell
          if(oldcell%n_points> 10000000) write (*,*) 'Split cell completed'
          call deletecell(oldcell)
    else
     cell => oldcell
     select case (cell%nsplitters)
     case (0)
       cell%cloudtype = 1
     case (1)
       cell%cloudtype = 2
     end select
     allocate(cell%loc(3,cell%n_points))
      allocate(cell%value(3,cell%n_points))
      cell%loc(:,1:cell%n_points) = current_cell_points_loc(:,1:cell%n_points)
      do n = 1, cell%n_points
        cell%value(ibase, n) = var_base(current_cell_points_loc(1,n), current_cell_points_loc(2,n), current_cell_points_loc(3,n))
        cell%value(itop, n) = var_top(current_cell_points_loc(1,n), current_cell_points_loc(2,n), current_cell_points_loc(3,n))
        cell%value(ivalue, n) = var_value(current_cell_points_loc(1,n), current_cell_points_loc(2,n), current_cell_points_loc(3,n))
        obj_mask(current_cell_points_loc(1,n), current_cell_points_loc(2,n), current_cell_points_loc(3,n)) = PROCESSED_OBJECT
      end do

    end if
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
    use constants, only: cbstep
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
    if (obj_mask(ii,j,t)==UNCLASSIFIED_IN_OBJECT .and. var_base(ii,j,t) < var_base(i,j,t) + cbstep) then
      obj_mask(ii,j,t) = from_point(4)
      num_boundary_points = num_boundary_points + 1
      boundary_points(1:3,num_boundary_points) = (/ii,j,t/)
      boundary_points(4,num_boundary_points) = from_point(4)
      num_points_per_new_cell(from_point(4)) = num_points_per_new_cell(from_point(4)) + 1
    end if
    !Look east
    ii = i + 1
    if (ii.gt.nx) ii = 1
    if (obj_mask(ii,j,t)==UNCLASSIFIED_IN_OBJECT .and. var_base(ii,j,t) < var_base(i,j,t) + cbstep) then
      obj_mask(ii,j,t) = from_point(4)
      num_boundary_points = num_boundary_points + 1
      boundary_points(1:3,num_boundary_points) = (/ii,j,t/)
      boundary_points(4,num_boundary_points) = from_point(4)
      num_points_per_new_cell(from_point(4)) = num_points_per_new_cell(from_point(4)) + 1
    end if
    !Look north
    jj = j - 1
    if (jj.le.0) jj = ny
    if (obj_mask(i,jj,t)==UNCLASSIFIED_IN_OBJECT .and. var_base(i,jj,t) < var_base(i,j,t) + cbstep) then
      obj_mask(i,jj,t) = from_point(4)
      num_boundary_points = num_boundary_points + 1
      boundary_points(1:3,num_boundary_points) = (/i,jj,t/)
      boundary_points(4,num_boundary_points) = from_point(4)
      num_points_per_new_cell(from_point(4)) = num_points_per_new_cell(from_point(4)) + 1
    end if
    !Look south
    jj = j + 1
    if (jj.gt.ny) jj = 1
    if (obj_mask(i,jj,t)==UNCLASSIFIED_IN_OBJECT .and. var_base(i,jj,t) < var_base(i,j,t) + cbstep) then
      obj_mask(i,jj,t) = from_point(4)
      num_boundary_points = num_boundary_points + 1
      boundary_points(1:3,num_boundary_points) = (/i,jj,t/)
      boundary_points(4,num_boundary_points) = from_point(4)
      num_points_per_new_cell(from_point(4)) = num_points_per_new_cell(from_point(4)) + 1
    end if
    !Look forward
    tt = t+1
    if (tt <= nt) then
      if (obj_mask(i,j,tt)==UNCLASSIFIED_IN_OBJECT .and. var_base(i,j,tt) < var_base(i,j,t) + cbstep) then
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
      if (obj_mask(i,j,tt)==UNCLASSIFIED_IN_OBJECT .and. var_base(i,j,tt) < var_base(i,j,t) + cbstep) then
        obj_mask(i,j,tt) = from_point(4)
        num_boundary_points = num_boundary_points + 1
        boundary_points(1:3,num_boundary_points) = (/i,j,tt/)
        boundary_points(4,num_boundary_points) = from_point(4)
        num_points_per_new_cell(from_point(4)) = num_points_per_new_cell(from_point(4)) + 1
      end if
    end if
  end subroutine findneighbour


  recursive subroutine newpassive(i, j, t, nr, totnewcells, nendlist, obj_mask, endlist, nractive)
    use tracking_common, only: tstart, nt, nx, ny
    use tracking_common, only: UNCLASSIFIED_IN_OBJECT

    integer(kind=2), intent(in) :: i, j, t
    integer, intent(in) :: totnewcells

    integer, intent(inout) :: nendlist
    integer, intent(inout), dimension(:,:)   :: endlist
    integer, intent(inout), dimension(:)     :: nr
    integer, intent(inout)     :: nractive
    integer(kind=4), intent(inout), dimension(:,:,:) :: obj_mask

    integer(kind=2) :: ii, jj, tt

    nr(totnewcells) = nr(totnewcells) + 1
!       if(nr(totnewcells) >= n_minparentel) nractive = 0
    nendlist = nendlist + 1
    obj_mask(i,j,t) = -1
    endlist(1:3,nendlist) = (/i,j,t/)
    endlist(4,nendlist)   = totnewcells

    !Look west
    ii = i - 1
    jj = j
    tt = t
    if (ii.le.0) ii = nx
    if (obj_mask(ii,jj,tt)==UNCLASSIFIED_IN_OBJECT) then
       call newpassive(i=ii, j=jj, t=tt, nr=nr, totnewcells=totnewcells, &
                       nendlist=nendlist, obj_mask=obj_mask, endlist=endlist, &
                       nractive=nractive)
    elseif (nractive > 0 .and. obj_mask(ii,jj,tt) /= nractive) then
      nractive = 0
    elseif (nractive <0 .and. obj_mask(ii,jj,tt) >0) then
      nractive = obj_mask(ii,jj,tt)
    end if
    !Look east
    ii = i + 1
    jj = j
    tt = t
    if (ii.gt.nx) ii = 1
    if (obj_mask(ii,jj,tt)==UNCLASSIFIED_IN_OBJECT) then
      call newpassive(i=ii,j=jj,t=tt, nr=nr, totnewcells=totnewcells, &
                      nendlist=nendlist, obj_mask=obj_mask, endlist=endlist, &
                      nractive=nractive)
    elseif (nractive > 0 .and. obj_mask(ii,jj,tt) /= nractive) then
      nractive = 0
    elseif (nractive <0 .and. obj_mask(ii,jj,tt) >0) then
      nractive = obj_mask(ii,jj,tt)
    end if

    !Look north
    ii = i
    jj = j - 1
    tt = t
    if (jj.le.0) jj = ny
    if (obj_mask(ii,jj,tt)==UNCLASSIFIED_IN_OBJECT) then
      call newpassive(i=ii, j=jj, t=tt, nr=nr, totnewcells=totnewcells, &
                      nendlist=nendlist, obj_mask=obj_mask, endlist=endlist, &
                      nractive=nractive)
    elseif (nractive > 0 .and. obj_mask(ii,jj,tt) /= nractive) then
      nractive = 0
    elseif (nractive <0 .and. obj_mask(ii,jj,tt) >0) then
      nractive = obj_mask(ii,jj,tt)
    end if

    !Look south
    ii = i
    jj = j + 1
    tt = t
    if (jj.gt.ny) jj = 1
    if (obj_mask(ii,jj,tt)==UNCLASSIFIED_IN_OBJECT) then
      call newpassive(i=ii, j=jj, t=tt, nr=nr, totnewcells=totnewcells, &
                      nendlist=nendlist, obj_mask=obj_mask, endlist=endlist, &
                      nractive=nractive)
    elseif (nractive > 0 .and. obj_mask(ii,jj,tt) /= nractive) then
      nractive = 0
    elseif (nractive <0 .and. obj_mask(ii,jj,tt) >0) then
      nractive = obj_mask(ii,jj,tt)
    end if

    !Look forward
    ii = i
    jj = j
    tt = t + 1
    if (tt <=nt) then
      if (obj_mask(ii,jj,tt)==UNCLASSIFIED_IN_OBJECT) then
        call newpassive(i=ii, j=jj, t=tt, nr=nr, totnewcells=totnewcells, &
                        nendlist=nendlist, obj_mask=obj_mask, endlist=endlist, &
                        nractive=nractive)
      elseif (nractive > 0 .and. obj_mask(ii,jj,tt) /= nractive) then
        nractive = 0
      elseif (nractive <0 .and. obj_mask(ii,jj,tt) >0) then
        nractive = obj_mask(ii,jj,tt)
      end if
    end if

    !Look backward
    ii = i
    jj = j
    tt = t - 1
    if (tt>=tstart) then
      if (obj_mask(ii,jj,tt)==UNCLASSIFIED_IN_OBJECT) then
        call newpassive(i=ii, j=jj, t=tt, nr=nr, totnewcells=totnewcells, &
                        nendlist=nendlist, obj_mask=obj_mask, endlist=endlist, &
                        nractive=nractive)
      elseif (nractive > 0 .and. obj_mask(ii,jj,tt) /= nractive) then
        nractive = 0
      elseif (nractive <0 .and. obj_mask(ii,jj,tt) >0) then
        nractive = obj_mask(ii,jj,tt)
      end if
    end if
  end subroutine newpassive
endmodule modtrack_cell_splitting
