module modtrack_cell_splitting
  use tracking_common, only: ibase, itop, ivalue
  use tracking_common, only: celltype, cellptr

  implicit none

  public splitcell

  private

  integer :: totcloudsystemnr = 0

  contains

  !> Find the parent cell for each element in the current cell and update
  !> `splitters` and `children` in each. Also creates a number of each
  !> parent cell which assigned to the last row of `list`. `nr` appears to
  !> count the number of "child" elements that are split by each cell.
  !> 
  !> the cloud mask at given point gets set to the number of splitters
  subroutine findnrsplitters_new(cell, parentarr, obj_mask, list, nr, nlist, current_cell_points_loc)
    use modarray, only: increase_array

    type(celltype),  pointer, intent(inout)                       :: cell
    type(cellptr),   allocatable, dimension(:,:,:), intent(inout) :: parentarr
    integer(kind=4), dimension(:,:,:), intent(inout)              :: obj_mask
    integer(kind=2), dimension(:,:), intent(in) :: current_cell_points_loc

    integer, dimension(:,:), intent(inout)   :: list
    integer, dimension(:), intent(inout)     :: nr
    integer, intent(inout)                   :: nlist

    integer :: nn, i, j, t, n
    logical :: lnewparent

    do nn = 1, cell%nelements
      !write(*,*) 'Find splitter for element',nn,'/',cell%nelements
      if (mod(nn,1000000) == 0) write(*,*) 'Find splitter for element',nn,'/',cell%nelements
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
        nlist = nlist + 1
        obj_mask(i,j,t) = n
        list(1,nlist) = i
        list(2,nlist) = j
        list(3,nlist) = t
        list(4,nlist) = n
        nr(n)         = nr(n) + 1
      end if
    end do
  end subroutine findnrsplitters_new
  

  !> Use a parent
  subroutine splitcell(cell, ncells, parentarr, obj_mask, var_base, var_top, var_value, current_cell_points_loc)
    use tracking_common, only: nrel_max
    use tracking_common, only: createcell
    use tracking_common, only: deletecell
    use tracking_common, only: print_cell_debug
    use tracking_common, only: MARKED_OBJECT, PROCESSED_OBJECT
    use constants, only: n_minparentel

    type(celltype), pointer, intent(inout)   :: cell
    integer, intent(inout)           :: ncells
    type(cellptr), allocatable, dimension(:,:,:), intent(inout) :: parentarr
    integer(kind=4), dimension(:,:,:), intent(inout) :: obj_mask

    integer(kind=2), dimension(:,:,:), intent(in) :: var_base
    integer(kind=2), dimension(:,:,:), intent(in) :: var_top
    integer(kind=2), dimension(:,:,:), intent(in) :: var_value
    integer(kind=2), dimension(:,:), intent(in) :: current_cell_points_loc

    integer, allocatable, dimension(:)     :: maxiter
    type(celltype), pointer                  :: oldcell
    type(cellptr), allocatable, dimension(:) :: newcells
    integer, allocatable, dimension(:,:)   :: list, newlist, endlist
    integer, allocatable, dimension(:)     :: nr
    integer :: n, tmin, tmax, iter
    integer :: nn, nlist, nnewlist, nendlist, npassive, totnewcells, nractive
    ! list has shape (4, cell%nelements) and stores (i,j,t,n) where i,j and the
    ! position indecies, t the time index and n the number of splitters for this
    ! element of the cell

    integer, allocatable, dimension(:,:)   :: split_neighbours_points_loc
    integer, allocatable, dimension(:)     :: split_neighbours_splitter_index
    integer, allocatable, dimension(:)     :: n_split_neighbours_points

    integer :: print_interval = 1

    !print *, "before splitcell"
    !call print_cell_debug(cell)

    oldcell => cell
    if(cell%nelements> 10000000) write (*,*) 'Begin Splitcell'

    allocate(split_neighbours_points_loc(3,oldcell%nelements))
    allocate(split_neighbours_splitter_index(oldcell%nelements))
    allocate(n_split_neighbours_points(oldcell%nelements))

    allocate(list(4,oldcell%nelements))
    allocate(newlist(4,oldcell%nelements))
    allocate(endlist(4,oldcell%nelements))
    allocate(nr(oldcell%nelements/2))

    nlist = 0
    nnewlist = 0
    newlist = 0
    nendlist = 0
    endlist = 0
    npassive = 0
    totnewcells = 0
    nr = 0

    tmin = minval(current_cell_points_loc(3,1:cell%nelements))
    tmax = maxval(current_cell_points_loc(3,1:cell%nelements))
    call findnrsplitters(cell=cell, parentarr=parentarr, obj_mask=obj_mask, list=list, &
                         nr=nr, nlist=nlist, current_cell_points_loc=current_cell_points_loc)
    if (oldcell%nsplitters >= 2) then
          if(oldcell%nelements> 10000000) then
            write (*,*) 'Split cell in ', cell%nsplitters, ' parts'
          endif

          allocate(maxiter(oldcell%nsplitters))

          ! insert `list` at current marker in `endlist`
          endlist(:,nendlist+1:nendlist+nlist) = list(:,1:nlist)
          nendlist = nendlist + nlist

          maxiter = max(nr(1:oldcell%nsplitters)/(20*n_minparentel),20)
          do iter = 1, maxval(maxiter)
            if (nlist == 0) exit
            if(mod(iter,100)==0) write (*,'(i12,a,i12,a,i12,a)') nlist, ' boundary points', &
                                                    oldcell%nelements - nendlist,'/' ,oldcell%nelements, ' points left'
            nnewlist = 0
            do nn =1, nlist
              !if (iter > maxiter(list(4,nn))) cycle
              call findneighbour(list=list, nn=nn, var_base=var_base, obj_mask=obj_mask, newlist=newlist, nr=nr, nnewlist=nnewlist)
            end do
            nlist = nnewlist
            list(:,1:nnewlist)  = newlist(:,1:nnewlist)
            endlist(:,nendlist+1:nendlist+nlist) = list(:,1:nlist)
            nendlist = nendlist + nlist
          end do
          nrel_max = max(nrel_max, 1)

          !Add cells for the number of parentcells
      !       npassive = 0
          totnewcells = oldcell%nsplitters
          if (nendlist < oldcell%nelements) then
            nlist = nendlist
            loop: do n = 1, oldcell%nelements
              if (obj_mask(current_cell_points_loc(1,n),&
                           current_cell_points_loc(2,n),&
                           current_cell_points_loc(3,n)&
                          ) == MARKED_OBJECT) then !Found new outflow region
      !             npassive = npassive + 1
                totnewcells = totnewcells + 1
                nr(totnewcells) = 0
                nractive = -1
                call newpassive(i=current_cell_points_loc(1,n), j=current_cell_points_loc(2,n), t=current_cell_points_loc(3,n), &
                                nr=nr, totnewcells=totnewcells, &
                                nendlist=nendlist, obj_mask=obj_mask, endlist=endlist, nractive=nractive)
                if (nractive > 0) then
                  endlist(4,1+nendlist-nr(totnewcells):nendlist) = nractive
                  nr(nractive) = nr(nractive) + nr(totnewcells)
                  nr(totnewcells) = 0
                  totnewcells = totnewcells -1
                end if
                if (mod(n,print_interval)==0) write (*,*) 'Passive cell',n,'/',oldcell%nelements
              end if
            end do loop
          end if

          allocate(newcells(totnewcells))
          do n = 1, oldcell%nsplitters
            if(mod(n,print_interval)==0) write (*,*) 'Create and fill new active cell', n, '/', oldcell%nsplitters

            call createcell(cell)
            newcells(n)%p => cell
            cell%nsplitters = 1
            allocate(cell%splitters(1))
            cell%splitters(1)%p => oldcell%splitters(n)%p
            cell%cloudtype = 4
            oldcell%splitters(n)%p%nchildren = 1
            deallocate(oldcell%splitters(n)%p%children)
            allocate(oldcell%splitters(n)%p%children(1))
            oldcell%splitters(n)%p%children(1)%p => cell
          end do
          do n = oldcell%nsplitters + 1, totnewcells
            if(mod(n,print_interval)==0) write (*,*) 'Create and fill new passive cell', n, '/', totnewcells
            call createcell(cell)
            newcells(n)%p => cell
            cell%cloudtype = 3
          end do
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
          do n = 1,totnewcells
           if (mod(n,print_interval)==0) write (*,*) 'Create split cell nr ', n,'/',totnewcells, 'with',nr(n),'elements'
      !         currid = oldcell%splitters(n)%p%id
            newcells(n)%p%nelements = nr(n)
            newcells(n)%p%cloudsystemnr = totcloudsystemnr
            allocate(newcells(n)%p%loc(3,nr(n)))
            allocate(newcells(n)%p%value(3,nr(n)))
          end do
          nr = 0
          do n = 1, nendlist
            if (mod(n,1000000)==0) write (*,*) 'Fill new cells ', n,'/',nendlist
            nn = endlist(4,n)
            nr(nn) = nr(nn) + 1
            newcells(nn)%p%loc(:,nr(nn)) = endlist(1:3,n)
            newcells(nn)%p%value(ibase,nr(nn)) = var_base(endlist(1,n),endlist(2,n),endlist(3,n))
            newcells(nn)%p%value(itop,nr(nn)) = var_top(endlist(1,n),endlist(2,n),endlist(3,n))
            newcells(nn)%p%value(ivalue,nr(nn)) = var_value(endlist(1,n),endlist(2,n),endlist(3,n))
          end do


          !Point at final new cell
          cell => newcells(totnewcells)%p
          do n = 1,totnewcells - 1
            ncells = ncells + 1
            newcells(n)%p%id = ncells
          end do
          !Remove the original cell
          if(oldcell%nelements> 10000000) write (*,*) 'Split cell completed'
          call deletecell(oldcell)
    else
     cell => oldcell
     select case (cell%nsplitters)
     case (0)
       cell%cloudtype = 1
     case (1)
       cell%cloudtype = 2
     end select
     allocate(cell%loc(3,cell%nelements))
      allocate(cell%value(3,cell%nelements))
      cell%loc(:,1:cell%nelements) = current_cell_points_loc(:,1:cell%nelements)
      do n = 1, cell%nelements
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
  subroutine findnrsplitters(cell, parentarr, obj_mask, list, nr, nlist, current_cell_points_loc)
    use modarray, only: increase_array

    type(celltype),  pointer, intent(inout)                       :: cell
    type(cellptr),   allocatable, dimension(:,:,:), intent(inout) :: parentarr
    integer(kind=4), dimension(:,:,:), intent(inout)              :: obj_mask
    integer(kind=2), dimension(:,:), intent(in) :: current_cell_points_loc

    integer, dimension(:,:), intent(inout)   :: list
    integer, dimension(:), intent(inout)     :: nr
    integer, intent(inout)                   :: nlist

    integer :: nn, i, j, t, n
    logical :: lnewparent

    do nn = 1, cell%nelements
      !write(*,*) 'Find splitter for element',nn,'/',cell%nelements
      if (mod(nn,1000000) == 0) write(*,*) 'Find splitter for element',nn,'/',cell%nelements
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
        nlist = nlist + 1
        obj_mask(i,j,t) = n
        list(1,nlist) = i
        list(2,nlist) = j
        list(3,nlist) = t
        list(4,nlist) = n
        nr(n)         = nr(n) + 1
      end if
    end do
  end subroutine findnrsplitters


  !> starting from position (in time and space) given as from (1,2,3) of the
  !> nnth element of `list` compare value at current position to neighbouring
  !> positions in time and space. The fourth component of `list` is used to
  !> assign into `obj_mask`
  subroutine findneighbour(list, nn, var_base, obj_mask, newlist, nr, nnewlist)
    use tracking_common, only: cbstep
    use tracking_common, only: tstart, nt, nx, ny
    use tracking_common, only: MARKED_OBJECT

    integer,         intent(in), dimension(:,:)      :: list
    integer,         intent(in)                      :: nn
    integer(kind=2), intent(in), dimension(:,:,:)    :: var_base

    integer(kind=4), intent(inout), dimension(:,:,:) :: obj_mask
    integer,         intent(inout), dimension(:,:)   :: newlist
    integer,         intent(inout), dimension(:)     :: nr
    integer,         intent(inout)                   :: nnewlist

    integer :: i, j, t
    integer :: ii ,jj, tt

    i = list(1,nn)
    j = list(2,nn)
    t = list(3,nn)

    !Look west
    ii = i - 1
    if (ii.le.0) ii = nx
    if (obj_mask(ii,j,t)==MARKED_OBJECT .and. var_base(ii,j,t) < var_base(i,j,t) + cbstep) then
      obj_mask(ii,j,t) = list(4,nn)
      nnewlist = nnewlist + 1
      newlist(1:3,nnewlist) = (/ii,j,t/)
      newlist(4,nnewlist) = list(4,nn)
      nr(list(4,nn)) = nr(list(4,nn)) + 1
    end if
    !Look east
    ii = i + 1
    if (ii.gt.nx) ii = 1
    if (obj_mask(ii,j,t)==MARKED_OBJECT .and. var_base(ii,j,t) < var_base(i,j,t) + cbstep) then
      obj_mask(ii,j,t) = list(4,nn)
      nnewlist = nnewlist + 1
      newlist(1:3,nnewlist) = (/ii,j,t/)
      newlist(4,nnewlist) = list(4,nn)
      nr(list(4,nn)) = nr(list(4,nn)) + 1
    end if
    !Look north
    jj = j - 1
    if (jj.le.0) jj = ny
    if (obj_mask(i,jj,t)==MARKED_OBJECT .and. var_base(i,jj,t) < var_base(i,j,t) + cbstep) then
      obj_mask(i,jj,t) = list(4,nn)
      nnewlist = nnewlist + 1
      newlist(1:3,nnewlist) = (/i,jj,t/)
      newlist(4,nnewlist) = list(4,nn)
      nr(list(4,nn)) = nr(list(4,nn)) + 1
    end if
    !Look south
    jj = j + 1
    if (jj.gt.ny) jj = 1
    if (obj_mask(i,jj,t)==MARKED_OBJECT .and. var_base(i,jj,t) < var_base(i,j,t) + cbstep) then
      obj_mask(i,jj,t) = list(4,nn)
      nnewlist = nnewlist + 1
      newlist(1:3,nnewlist) = (/i,jj,t/)
      newlist(4,nnewlist) = list(4,nn)
      nr(list(4,nn)) = nr(list(4,nn)) + 1
    end if
    !Look forward
    tt = t+1
    if (tt <= nt) then
      if (obj_mask(i,j,tt)==MARKED_OBJECT .and. var_base(i,j,tt) < var_base(i,j,t) + cbstep) then
        obj_mask(i,j,tt) = list(4,nn)
        nnewlist = nnewlist + 1
        newlist(1:3,nnewlist) = (/i,j,tt/)
        newlist(4,nnewlist) = list(4,nn)
        nr(list(4,nn)) = nr(list(4,nn)) + 1
      end if
    end if
    !Look backward
    tt = t - 1
    if (tt >=tstart) then
      if (obj_mask(i,j,tt)==MARKED_OBJECT .and. var_base(i,j,tt) < var_base(i,j,t) + cbstep) then
        obj_mask(i,j,tt) = list(4,nn)
        nnewlist = nnewlist + 1
        newlist(1:3,nnewlist) = (/i,j,tt/)
        newlist(4,nnewlist) = list(4,nn)
        nr(list(4,nn)) = nr(list(4,nn)) + 1
      end if
    end if
  end subroutine findneighbour


  recursive subroutine newpassive(i, j, t, nr, totnewcells, nendlist, obj_mask, endlist, nractive)
    use tracking_common, only: tstart, nt, nx, ny
    use tracking_common, only: MARKED_OBJECT

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
    if (obj_mask(ii,jj,tt)==MARKED_OBJECT) then
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
    if (obj_mask(ii,jj,tt)==MARKED_OBJECT) then
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
    if (obj_mask(ii,jj,tt)==MARKED_OBJECT) then
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
    if (obj_mask(ii,jj,tt)==MARKED_OBJECT) then
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
      if (obj_mask(ii,jj,tt)==MARKED_OBJECT) then
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
      if (obj_mask(ii,jj,tt)==MARKED_OBJECT) then
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
