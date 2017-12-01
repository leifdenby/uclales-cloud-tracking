module modtrack
  use tracking_common, only: celltype, cellptr

  use tracking_common, only: dt, dx, dy
  use tracking_common, only: ibase, itop, ivalue
  use tracking_common, only: nrel_max, nt, nx, ny
  use tracking_common, only: tstart
  use tracking_common, only: minparentel
  use tracking_common, only: cbstep

  use tracking_common, only: bool
  use tracking_common, only: var

  use modarray, only: increase_array

  implicit none

  private

  integer :: number
  integer :: totcloudsystemnr = 0
  integer(kind=2), dimension(:,:), allocatable :: cellloc

  ! main
  public dotracking, fillparentarr, findparents

  ! modstatistics
  public nextcell, firstcell

  contains

  !> Given the input `cell` create a new `celltype` instance with the same
  !>`head` and `next`, with `previous` pointing to the current cell and set
  !>`next` of the current `cell` to the new `celltype` instance.
  !
  !! @ returns (through reassignment) a pointer to the new `celltype` instance
  subroutine createcell(cell)
    type(celltype), pointer, intent(inout) :: cell
    type(celltype), pointer       :: tmp

    if (associated(cell)) then
      allocate(tmp)
      tmp%previous => cell
      tmp%next     => cell%next
      tmp%head     => cell%head
      cell%next    => tmp
      cell => tmp
    else
      allocate(cell)
      cell%head => cell
      cell%previous => null()
      cell%next     => null()
    end if

    cell%nelements = 0
    cell%nchildren = 0
    cell%nsplitters= 0
    cell%nparents  = 0
    cell%cloudsystemnr = 0
    cell%cloudtype = -1
  end subroutine createcell

  !> Given `cell` remove all it and all associated allocated "splitters",
  !> "parents", "children"
  !
  !! The `next`, `previous` linked-list pointers are re-assigned so that a new
  !! contigous link-list is created, and a pointer into that list (either the
  !! "next" or "previous" element) is re-assigned to the `cell` argument
  subroutine deletecell(cell)
    type(celltype), pointer, intent(inout) :: cell
    type(celltype), pointer       :: tmp
    integer :: ierr
    tmp => NULL()
    if (associated(cell%previous)) then
      tmp  => cell%previous
      tmp%next=> cell%next
    elseif (associated(cell%next)) then
      tmp  => cell%next
      tmp%previous=> cell%previous
    end if
    deallocate(cell%loc, cell%value, cell%splitters, cell%parents, cell%children, stat=ierr)
    deallocate(cell)
    if (associated(tmp)) then
      cell=>tmp
    end if
  end subroutine deletecell

  !> Through using `deletecell` de-allocate all cells (and associated "parents",
  !> etc) in the `next`, `previous` linked-list that the argument `cell` exists in
  subroutine delete_all(cell)
    type(celltype), pointer, intent(inout) :: cell
    integer :: iret
    iret = firstcell(cell)
    iret = 0
    do
    iret = iret + 1
      if (.not. associated(cell)) exit
      call deletecell(cell)
    end do
  end subroutine delete_all

  !> Assign the `cell` argument to the "next" `celltype` attribute of the supplied cell if the
  !> supplied cell is pointing to something
  !! @returns 0: `cell` is associated, -1: otherwise
  integer function nextcell(cell)
    type(celltype), pointer, intent(inout) :: cell

    nextcell = -1
    if (associated(cell%next)) then
      cell => cell%next
      nextcell = 0
    end if

  end function nextcell

  !> Assign the `cell` argument to the "head" `celltype` attribute of the supplied cell if the
  !supplied cell is pointing to something
  !@ return 0: `cell` is associated, -1: otherwise
  integer function firstcell(cell)
    type(celltype), pointer, intent(inout) :: cell
    firstcell = -1
    if(associated(cell)) then
      cell => cell%head
      firstcell = 0
    end if
  end function firstcell

  subroutine splitcell(cell, ncells, parentarr)
    type(celltype), pointer, intent(inout)   :: cell
    integer, intent(inout)           :: ncells
    type(cellptr), allocatable, dimension(:,:,:), intent(inout) :: parentarr
    type(celltype), pointer                  :: oldcell
    type(cellptr), allocatable, dimension(:) :: newcells
    logical :: lnewparent
    integer, allocatable, dimension(:,:)   :: list, newlist, endlist
    integer, allocatable, dimension(:)     :: nr, maxiter
    integer :: i, j, t, n, nnn, ii, jj, tt, currid, tmin, tmax, iter, maxpassive=0
    integer :: nn, nlist, nnewlist, nendlist, npassive, totnewcells, nractive
    ! list has shape (4, cell%nelements) and stores (i,j,t,n) where i,j and the
    ! position indecies, t the time index and n the number of splitters for this
    ! element of the cell

    oldcell => cell
    if(cell%nelements> 10000000) write (*,*) 'Begin Splitcell'

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

    tmin = minval(cellloc(3,1:cell%nelements))
    tmax = maxval(cellloc(3,1:cell%nelements))
    call findnrsplitters
    if (oldcell%nsplitters >= 2) then
      if(oldcell%nelements> 10000000) write (*,*) 'Split cell in ', cell%nsplitters, ' parts'
      allocate(maxiter(oldcell%nsplitters))

      endlist(:,nendlist+1:nendlist+nlist) = list(:,1:nlist)
      nendlist = nendlist + nlist

      maxiter = max(nr(1:oldcell%nsplitters)/(20*minparentel),5)
      do iter = 1, maxval(maxiter)
        if (nlist == 0) exit
        if(mod(iter,100)==0) write (*,'(i12,a,i12,a,i12,a)') nlist, ' boundary points', &
                                                oldcell%nelements - nendlist,'/' ,oldcell%nelements, ' points left'
        nnewlist = 0
        do nn =1, nlist
          if (iter > maxiter(list(4,nn))) cycle
          call findneighbour
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
          if (bool(cellloc(1,n), cellloc(2,n), cellloc(3,n)) == -2) then !Found new outflow region
!             npassive = npassive + 1
            totnewcells = totnewcells + 1
            nr(totnewcells) = 0
            nractive = -1
            call newpassive(cellloc(1,n), cellloc(2,n), cellloc(3,n))
            if (nractive > 0) then
              endlist(4,1+nendlist-nr(totnewcells):nendlist) = nractive
              nr(nractive) = nr(nractive) + nr(totnewcells)
              nr(totnewcells) = 0
              totnewcells = totnewcells -1
            end if
            if (mod(n,1000)==0) write (*,*) 'Passive cell',n,'/',oldcell%nelements
          end if
        end do loop
      end if

      allocate(newcells(totnewcells))
      do n = 1, oldcell%nsplitters
        if(mod(n,1000)==0) write (*,*) 'Create and fill new active cell', n, '/', oldcell%nsplitters

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
        if(mod(n,1000)==0) write (*,*) 'Create and fill new passive cell', n, '/', totnewcells
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
       if (mod(n,1000)==0) write (*,*) 'Create split cell nr ', n,'/',totnewcells, 'with',nr(n),'elements'
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
        newcells(nn)%p%value(1,nr(nn)) = var(endlist(1,n),endlist(2,n),endlist(3,n),ibase)
        newcells(nn)%p%value(2,nr(nn)) = var(endlist(1,n),endlist(2,n),endlist(3,n),itop)
        newcells(nn)%p%value(3,nr(nn)) = var(endlist(1,n),endlist(2,n),endlist(3,n),ivalue)
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
      cell%loc(:,1:cell%nelements) = cellloc(:,1:cell%nelements)
      do n = 1, cell%nelements
        cell%value(1, n) = var(cellloc(1,n), cellloc(2,n), cellloc(3,n), ibase)
        cell%value(2, n) = var(cellloc(1,n), cellloc(2,n), cellloc(3,n), itop)
        cell%value(3, n) = var(cellloc(1,n), cellloc(2,n), cellloc(3,n), ivalue)
        bool(cellloc(1,n), cellloc(2,n), cellloc(3,n)) = -1
      end do

    end if
  contains
    !> Find the parent cell for each element in the current cell and update
    !> `splitters` and `children` in each. Also creates a number of each
    !> parent cell which assigned to the last row of `list`. `nr` appears to
    !> count the number of "child" elements that are split by each cell.
    subroutine findnrsplitters
      do nn = 1, cell%nelements
  !      write(*,*) 'Find splitter for element',nn,'/',cell%nelements
        if (mod(nn,1000000) == 0) write(*,*) 'Find splitter for element',nn,'/',cell%nelements
        i = cellloc(1,nn)
        j = cellloc(2,nn)
        t = cellloc(3,nn)
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
          bool(i,j,t) = n
          list(1,nlist) = i
          list(2,nlist) = j
          list(3,nlist) = t
          list(4,nlist) = n
          nr(n)         = nr(n) + 1
        end if
      end do
    end subroutine findnrsplitters
    subroutine findneighbour
      i = list(1,nn)
      j = list(2,nn)
      t = list(3,nn)

      !Look west
      ii = i - 1
      if (ii.le.0) ii = nx
      if (bool(ii,j,t)==-2 .and. var(ii,j,t,ibase) < var(i,j,t,ibase) + cbstep) then
        bool(ii,j,t) = list(4,nn)
        nnewlist = nnewlist + 1
        newlist(1:3,nnewlist) = (/ii,j,t/)
        newlist(4,nnewlist) = list(4,nn)
        nr(list(4,nn)) = nr(list(4,nn)) + 1
      end if
      !Look east
      ii = i + 1
      if (ii.gt.nx) ii = 1
      if (bool(ii,j,t)==-2 .and. var(ii,j,t,ibase) < var(i,j,t,ibase) + cbstep) then
        bool(ii,j,t) = list(4,nn)
        nnewlist = nnewlist + 1
        newlist(1:3,nnewlist) = (/ii,j,t/)
        newlist(4,nnewlist) = list(4,nn)
        nr(list(4,nn)) = nr(list(4,nn)) + 1
      end if
      !Look north
      jj = j - 1
      if (jj.le.0) jj = ny
      if (bool(i,jj,t)==-2 .and. var(i,jj,t,ibase) < var(i,j,t,ibase) + cbstep) then
        bool(i,jj,t) = list(4,nn)
        nnewlist = nnewlist + 1
        newlist(1:3,nnewlist) = (/i,jj,t/)
        newlist(4,nnewlist) = list(4,nn)
        nr(list(4,nn)) = nr(list(4,nn)) + 1
      end if
      !Look south
      jj = j + 1
      if (jj.gt.ny) jj = 1
      if (bool(i,jj,t)==-2 .and. var(i,jj,t,ibase) < var(i,j,t,ibase) + cbstep) then
        bool(i,jj,t) = list(4,nn)
        nnewlist = nnewlist + 1
        newlist(1:3,nnewlist) = (/i,jj,t/)
        newlist(4,nnewlist) = list(4,nn)
        nr(list(4,nn)) = nr(list(4,nn)) + 1
      end if
      !Look forward
      tt = t+1
      if (tt <= nt) then
        if (bool(i,j,tt)==-2 .and. var(i,j,tt,ibase) < var(i,j,t,ibase) + cbstep) then
          bool(i,j,tt) = list(4,nn)
          nnewlist = nnewlist + 1
          newlist(1:3,nnewlist) = (/i,j,tt/)
          newlist(4,nnewlist) = list(4,nn)
          nr(list(4,nn)) = nr(list(4,nn)) + 1
        end if
      end if
      !Look backward
      tt = t - 1
      if (tt >=tstart) then
        if (bool(i,j,tt)==-2 .and. var(i,j,tt,ibase) < var(i,j,t,ibase) + cbstep) then
          bool(i,j,tt) = list(4,nn)
          nnewlist = nnewlist + 1
          newlist(1:3,nnewlist) = (/i,j,tt/)
          newlist(4,nnewlist) = list(4,nn)
          nr(list(4,nn)) = nr(list(4,nn)) + 1
        end if
      end if
    end subroutine findneighbour
    recursive subroutine newpassive(i, j, t)
      integer(kind=2), intent(in) :: i, j, t
      integer(kind=2) :: ii, jj, tt


      nr(totnewcells) = nr(totnewcells) + 1
!       if(nr(totnewcells) >= minparentel) nractive = 0
      nendlist = nendlist + 1
      bool(i,j,t) = -1
      endlist(1:3,nendlist) = (/i,j,t/)
      endlist(4,nendlist)   = totnewcells

      !Look west
      ii = i - 1
      jj = j
      tt = t
      if (ii.le.0) ii = nx
      if (bool(ii,jj,tt)==-2) then
        call newpassive(ii,jj,tt)
      elseif (nractive > 0 .and. bool(ii,jj,tt) /= nractive) then
        nractive = 0
      elseif (nractive <0 .and. bool(ii,jj,tt) >0) then
        nractive = bool(ii,jj,tt)
      end if
      !Look east
      ii = i + 1
      jj = j
      tt = t
      if (ii.gt.nx) ii = 1
      if (bool(ii,jj,tt)==-2) then
        call newpassive(ii,jj,tt)
      elseif (nractive > 0 .and. bool(ii,jj,tt) /= nractive) then
        nractive = 0
      elseif (nractive <0 .and. bool(ii,jj,tt) >0) then
        nractive = bool(ii,jj,tt)
      end if

      !Look north
      ii = i
      jj = j - 1
      tt = t
      if (jj.le.0) jj = ny
      if (bool(ii,jj,tt)==-2) then
        call newpassive(ii,jj,tt)
      elseif (nractive > 0 .and. bool(ii,jj,tt) /= nractive) then
        nractive = 0
      elseif (nractive <0 .and. bool(ii,jj,tt) >0) then
        nractive = bool(ii,jj,tt)
      end if

      !Look south
      ii = i
      jj = j + 1
      tt = t
      if (jj.gt.ny) jj = 1
      if (bool(ii,jj,tt)==-2) then
        call newpassive(ii,jj,tt)
      elseif (nractive > 0 .and. bool(ii,jj,tt) /= nractive) then
        nractive = 0
      elseif (nractive <0 .and. bool(ii,jj,tt) >0) then
        nractive = bool(ii,jj,tt)
      end if

      !Look forward
      ii = i
      jj = j
      tt = t + 1
      if (tt <=nt) then
        if (bool(ii,jj,tt)==-2) then
          call newpassive(ii,jj,tt)
        elseif (nractive > 0 .and. bool(ii,jj,tt) /= nractive) then
          nractive = 0
        elseif (nractive <0 .and. bool(ii,jj,tt) >0) then
          nractive = bool(ii,jj,tt)
        end if
      end if

      !Look backward
      ii = i
      jj = j
      tt = t - 1
      if (tt>=tstart) then
        if (bool(ii,jj,tt)==-2) then
          call newpassive(ii,jj,tt)
        elseif (nractive > 0 .and. bool(ii,jj,tt) /= nractive) then
          nractive = 0
        elseif (nractive <0 .and. bool(ii,jj,tt) >0) then
          nractive = bool(ii,jj,tt)
        end if
      end if
    end subroutine newpassive
  end subroutine splitcell


  subroutine findparents(cell, parentarr, base, top)
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
      do nn = 1, cell%nelements
        i = cell%loc(1,nn)
        j = cell%loc(2,nn)
        t = cell%loc(3,nn)
        if (associated(parentarr(i,j,t)%p)) then
          if (base(i,j,t) <= cell%value(itop,nn) .and. &
              top(i,j,t)  >= cell%value(ibase,nn) ) then
            lnewparent = .true.
            do n = 1, cell%nparents
              if (cell%parents(n)%p%id == parentarr(i,j,t)%p%id) then
                lnewparent = .false.
                exit
              end if
            end do
            if (lnewparent) then
              cell%nparents = cell%nparents + 1
              nrel_max = max(nrel_max, cell%nparents)
              call increase_array(cell%parents, cell%nparents)
              cell%parents(cell%nparents)%p => parentarr(i,j,t)%p

              parentarr(i,j,t)%p%nchildren = parentarr(i,j,t)%p%nchildren + 1
              nrel_max = max(nrel_max, parentarr(i,j,t)%p%nchildren)
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

  !> Copy the information from the temporary array `cellloc` into the `cell`
  !> instance's `loc` (location) and `value` attributes.
  !!
  !! - `cellloc` contains all the positions in space of time of datapoints which
  !! are part of the current cell
  !! - copy into `cell.loc` the location in space of time of each element
  !! - set in `cell.value` the scalar values of each element (which a defined in
  !! both space and time) wrt cloud-top, cloud-base and "ivalue" ??
  !! - bool is set to `-1` for all datapoints for the elements in the current
  !! cell
  !!
  !! @TODO what does `ivalue` mean here?
  subroutine finalizecell(cell, ncells, parentarr)
    type(celltype), pointer, intent(inout)                   :: cell
    integer, intent(out)                             :: ncells
    type(cellptr), allocatable, dimension(:,:,:), intent(inout), optional :: parentarr
    integer :: n, i, j, t, nr

    if (present(parentarr)) then
      call splitcell(cell, ncells, parentarr)
    else
      allocate(cell%loc(3,cell%nelements))
      allocate(cell%value(3,cell%nelements))
      cell%loc(:,1:cell%nelements) = cellloc(:,1:cell%nelements)
      do n = 1, cell%nelements
        cell%value(1, n) = var(cellloc(1,n), cellloc(2,n), cellloc(3,n), ibase)
        cell%value(2, n) = var(cellloc(1,n), cellloc(2,n), cellloc(3,n), itop)
        cell%value(3, n) = var(cellloc(1,n), cellloc(2,n), cellloc(3,n), ivalue)
        bool(cellloc(1,n), cellloc(2,n), cellloc(3,n)) = -1
      end do
    end if
    ncells = ncells + 1
    cell%id = ncells
  end subroutine finalizecell

  !> Given the current location (i,j,t) in space and time look at neighbouring
  !> points in space and time and if they satisfy the constraints from being part of
  !> the same cell:
  !!
  !! set bool=-2 so that this data-point is not considered twice for multiple
  !! cells, store the position in space and time into the `cellloc` array and
  !! increment the `nelements` counter on the provided cell
  !!
  !! TODO: Looking west/east/north/south and truncating the indexing near the
  !! edge is unnecessarily costly, the current element will be checked twice
  recursive subroutine newelement(i, j, t, cell)
    integer, intent(in)                                      :: i, j, t
    type(celltype),pointer, intent(inout)                    :: cell
    integer :: ii, jj, tt

    cell%nelements = cell%nelements + 1
    if (mod(cell%nelements,10000000) == 0) write(*,*) 'Cell element ', cell%nelements
    cellloc(1,cell%nelements) = i
    cellloc(2,cell%nelements) = j
    cellloc(3,cell%nelements) = t

    bool(i, j, t) = -2

    !Look west
    ii = i - 1
    jj = j
    tt = t
    if (ii.le.0) ii = nx
    if (bool(ii,jj,tt)==0) then
      if (var(i,j,t,ibase) <= var(ii,jj,tt,itop) .and. var(i,j,t,itop) >= var(ii,jj,tt,ibase)) then
        call newelement(ii,jj,tt,cell)
      end if
    end if
    !Look east
    ii = i + 1
    jj = j
    tt = t
    if (ii.gt.nx) ii = 1
    if (bool(ii,jj,tt)==0) then
      if (var(i,j,t,ibase) <= var(ii,jj,tt,itop) .and. var(i,j,t,itop) >= var(ii,jj,tt,ibase)) then
        call newelement(ii,jj,tt,cell)
      end if
    end if

    !Look north
    ii = i
    jj = j - 1
    tt = t
    if (jj.le.0) jj = ny
    if (bool(ii,jj,tt)==0) then
      if (var(i,j,t,ibase) <= var(ii,jj,tt,itop) .and. var(i,j,t,itop) >= var(ii,jj,tt,ibase)) then
        call newelement(ii,jj,tt,cell)
      end if
    end if

    !Look south
    ii = i
    jj = j + 1
    tt = t
    if (jj.gt.ny) jj = 1
    if (bool(ii,jj,tt)==0) then
      if (var(i,j,t,ibase) <= var(ii,jj,tt,itop) .and. var(i,j,t,itop) >= var(ii,jj,tt,ibase)) then
        call newelement(ii,jj,tt,cell)
      end if
    end if

    !Look forward
    ii = i
    jj = j
    tt =t + 1
    if (tt <= nt) then
      if (bool(ii,jj,tt)==0) then
        if (var(i,j,t,ibase) <= var(ii,jj,tt,itop) .and. var(i,j,t,itop) >= var(ii,jj,tt,ibase)) then
          call newelement(ii,jj,tt,cell)
        end if
      end if
    end if
    !Look backward
    ii = i
    jj = j
    tt = t - 1
    if (tt >= tstart) then
      if (bool(ii,jj,tt)==0) then
        if (var(i,j,t,ibase) <= var(ii,jj,tt,itop) .and. var(i,j,t,itop) >= var(ii,jj,tt,ibase)) then
          call newelement(ii,jj,tt,cell)
        end if
      end if
    end if
  end subroutine newelement

  !> Iterate over all data-points in space and time and construct `celltype`
  !> instances (stored in a linked list through the `next`/`previous` attributes) 
  subroutine dotracking(cell, ncells, nmincells, parentarr)
    type(celltype), pointer, intent(inout)                   :: cell
    integer, intent(out)                             :: ncells
    integer, intent(in)                                      :: nmincells
    type(cellptr), allocatable, dimension(:,:,:), intent(inout), optional :: parentarr
    integer :: i, j, t
    write (*,*) '.. entering tracking'

    allocate(cellloc(3,ceiling(min(0.3*(huge(1)-2),0.5*real(nx)*real(ny)*real(nt-tstart)))))
print *, 'cellloc',shape(cellloc),0.3*huge(1), 0.5*real(nx)*real(ny)*real(nt-tstart)
    nullify(cell)
    do t = tstart, nt
      if(mod(t,10)==0) write (*,'(A,I10,A,I10)') "Time =", t,"  Nr of Cells", ncells
      do  j = 1, ny
        do i = 1, nx
          if (bool(i,j,t)==0) then
!             write (*,*) 'Create a new cell'
            call createcell(cell)
            number = 0
            call newelement(i, j, t, cell)
            if (cell%nelements >= nmincells) then
              if(cell%nelements> 10000000) write (*,*) '..finalizing cell'
              call finalizecell(cell, ncells, parentarr)
            else
              call deletecell(cell)
            end if
          end if
        end do
      end do
    end do
    write (*,*) '.. leaving tracking'
    deallocate(cellloc)
  end subroutine dotracking

  !> Create a mapping from all datapoints in space and time to the cell
  !>  associated with each datapoint, but only for cells where the minimum "cloud-base"
  !>  value is less than minimum base height provided (this is actually and array
  !>  that spans all timesteps, but only the first time-step index of first
  !>  element in the cell is used).
  !>
  !! Also optionally store the "cloud-base" and "cloud-top" value into `base` and
  !! `top` arrays
  !!
  !! @TODO Does this subroutine actually read from `parentarr`?
  subroutine fillparentarr(cell, minbase, parentarr, base, top)
    use modnetcdf, only : fillvalue_i16
    type(celltype), pointer, intent(inout)         :: cell
    type(cellptr), allocatable, dimension(:,:,:), intent(inout)   :: parentarr
    integer(kind=2), allocatable,dimension(:), intent(in)      :: minbase
    integer(kind=2), allocatable,dimension(:,:,:), intent(inout), optional  :: base, top
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
      if (cell%nelements > minparentel) then
        if (minval(cell%value(ibase,1:cell%nelements))< minbase(cell%loc(3,1))) then
          do n = 1, cell%nelements
            parentarr(cell%loc(1,n),cell%loc(2,n),cell%loc(3,n))%p => cell
            if (present(base)) then
              base(cell%loc(1,n),cell%loc(2,n),cell%loc(3,n))        =  cell%value(ibase,n)
              top(cell%loc(1,n),cell%loc(2,n),cell%loc(3,n))        =  cell%value(itop,n)
            end if
          end do
        end if
      end if
      iret =  nextcell(cell)
    end do
    write (*,*) '.. leaving fillparentarr'
  end subroutine fillparentarr


  !> Check if successive timesteps appear to strange number of `bool==0`
  !> data-points
  !!
  !! Strange defined as either:
  !! - number of cells with bool==0 was above 100 in previous step and is 0 in current 
  !! - succesive timesteps has less than 25% of bool==0 than previous timestep
  !!
  !! There's something strange whereby bool is actually changed if either of the
  !! above is true, I think this is to always relate to the last timestep that
  !! was deemed to be good for comparison, maybe...
  subroutine checkframes
    integer :: t, count1, count2
    real    :: hlp
      !..check for empty frames
    write(0,*) 'Check frames'
    count2 = 0
    do t = 2, nt
      count1 = count(bool(:,:,t)==0)
      hlp = float(count1)/float(count2)
      if (count1.eq.0.and.count2.gt.100) then
          write (0,*) 'WARNING: empty frame at t = ',t
          bool(:,:,t)= bool(:,:,t-1)
          count1 = count2
      elseif (hlp.lt.0.25) then
          write (0,*) 'WARNING: strange frame at t = ',t
          bool(:,:,t)= bool(:,:,t-1)
          count1 = count2
      end if
      count2 = count1
    end do

  end subroutine checkframes
end module modtrack
