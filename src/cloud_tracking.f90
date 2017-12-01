
module modtrack
  implicit none
  !Declare

  type :: cellptr
    type(celltype), pointer :: p
  end type cellptr

  type :: celltype
    integer :: id
    integer :: nelements
    integer(kind=2), allocatable, dimension(:,:)     :: loc
    integer(kind=2), allocatable, dimension (:,:)       :: value
    integer                                  :: nsplitters
    type(cellptr), allocatable, dimension(:) :: splitters
    integer                                  :: nparents
    type(cellptr), allocatable, dimension(:) :: parents
    integer                                  :: nchildren
    type(cellptr), allocatable, dimension(:) :: children
    integer                                  :: cloudsystemnr
    integer                                  :: cloudtype
!     type(cellptr), allocatable, dimension(:) :: siblings
    type(celltype), pointer :: previous
    type(celltype), pointer :: next
    type(celltype), pointer :: head
  end type celltype


  type(celltype), pointer :: core, cloud, rain, thermal
  integer :: ncores, nclouds, nrains, nthermals, nrel_max
  integer :: nx, ny, nt, nvar, tstart
  integer(kind=2), dimension(:,:,:,:), allocatable :: var
  integer(kind=4), dimension(:,:,:), allocatable :: bool
  integer(kind=2), dimension(:,:), allocatable :: cellloc
  integer :: number
  integer :: totcloudsystemnr = 0
  integer, parameter :: ibase = 1, itop = 2
  integer :: ivalue = 3, icore = -1, ilwp = -1, irain = -1, ithermal = -1
  real    :: dx, dy, dt
  integer :: minparentel
  integer(kind=2)    :: cbstep
  logical, parameter :: ldebug = .false., lsiblings = .false.
  interface increase_array
    module procedure increase_array_i
    module procedure increase_array_r
    module procedure increase_array_p
  end interface increase_array

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

  !> Increase an array of reals by allocating a new array of size `newsize` and copy the old array to the beginning of the new if
  !> the old one was allocated
  subroutine increase_array_r(array, newsize)
    real,    dimension(:,:), allocatable, intent(inout) :: array
    integer, dimension(:), intent(in)      :: newsize
    real,    dimension(:,:), allocatable   :: tmp

    if(allocated(array)) call move_alloc(array, tmp)
    allocate (array(newsize(1), newsize(2)))
    if (allocated(tmp)) then
      array(1:size(tmp,1),1:size(tmp,2)) = tmp
    end if

  end subroutine increase_array_r

  !> Increase an array of integers by allocating a new array of size `newsize` and copy the old array to the beginning of the new if
  !> the old one was allocated
  subroutine increase_array_i(array, newsize)
    integer, dimension(:,:),  allocatable, intent(inout) :: array
    integer, dimension(:), intent(in)      :: newsize
    integer, dimension(:,:), allocatable   :: tmp

    if(allocated(array)) call move_alloc(array, tmp)
    allocate (array(newsize(1), newsize(2)))
    if (allocated(tmp)) then
      array(1:size(tmp,1),1:size(tmp,2)) = tmp
    end if

  end subroutine increase_array_i

  !> Increase an array of the `cellptr` datatype by allocating a new array of size `newsize` and copy the old array to the beginning of the new if
  !> the old one was allocated
  subroutine increase_array_p(array, newsize)
    type(cellptr), dimension(:),  allocatable, intent(inout) :: array
    integer, intent(in)                        :: newsize
    type(cellptr), dimension(:),  allocatable  :: tmp

    if(allocated(array)) call move_alloc(array, tmp)
    allocate (array(newsize))
    if (allocated(tmp)) then
      array(1:size(tmp)) = tmp
    end if

  end subroutine increase_array_p

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

  subroutine dostatistics(cell, ncells, icell, fid, ivar, heightzero, heightrange, maxheight, valzero, valrange, ovarstem)
    use modnetcdf
    type(celltype), pointer, intent(inout)    :: cell
    integer, intent(in)               :: ncells
    integer, intent(in)                       :: fid, icell
    real, intent(in)                          :: heightzero, heightrange, maxheight,valzero, valrange
    type(netcdfvar), dimension(:), intent(in) :: ivar
    type(netcdfvar)              , intent(in) :: ovarstem
    type(netcdfvar) :: ovar, nrnc, timenc, relnc
    integer :: i, j, t, tt,kk, n, nn,  tmin, tmax, cummt, iret, idev, jdev
    integer :: nel
    integer, parameter :: nbuckets = 3
    real, dimension(nbuckets) :: r_bucketsize
    integer, dimension(nbuckets) :: bucketsize, bucket_min, bucket_max
    character(80), dimension(nbuckets) :: bucketname, bucketlongname
    integer, allocatable, dimension(:) :: tlength, tdistr
    integer, allocatable, dimension(:) ::  icenter, jcenter, ianchor, janchor, npts
    integer, dimension(nx, ny, tstart:nt) :: slab
    real, dimension(:,:), allocatable :: base, top, area, vol, val, xcenter, ycenter, maxarea, maxarealoc, recon
    real, dimension(:), allocatable   :: duration, mintime, maxtime
    real :: dz = 25, rbase, rtop
    integer, dimension(:,:), allocatable :: relatives
    integer, dimension(:), allocatable   :: id, nrelatives
real :: time
    write (*,*) '.. entering statistics'

    !Loop over the cells  - fill the cell-length distribution and the xyt slab
    slab(1:nx, 1:ny, tstart:nt) = 0
    slab(1:nx, 1:ny, tstart:nt) = -1000000000
    slab(1:nx, 1:ny, tstart:nt) = fillvalue_i
    bucketsize = 0
    bucket_max = 0
    bucket_min = 0
    bucketname = (/'sm','la','hu'/)
    bucketlongname = (/'Small','Large','Huge '/)

    allocate(nrnc%dim(1))
    allocate(nrnc%dimids(1))
    allocate(timenc%dim(1))
    allocate(timenc%dimids(1))
    allocate(relnc%dim(1))
    allocate(relnc%dimids(1))


    !Write nr of relatives dimension
    relnc%name     = trim(ovarstem%name)//'rel'
    relnc%longname = trim(ovarstem%name)// ' relatives'
    relnc%units    = '#'
    relnc%dim(1)   = nrel_max
    call define_ncdim(fid, relnc, nf90_int)
    call write_ncvar(fid, relnc, (/(i, i=1,relnc%dim(1))/))

    n = 0
    allocate(tlength(ncells))
    iret = firstcell(cell)
    do
      if (iret == -1) exit
      n = n + 1
      tmin = minval(cell%loc(3,1:cell%nelements))
      tmax = maxval(cell%loc(3,1:cell%nelements))
      tlength(n) = tmax - tmin + 1
      iret = nextcell(cell)
    end do
    allocate(tdistr(maxval(tlength)))
    tdistr = 0
    cummt  = 0
    r_bucketsize = (/1.,1., 1./) * ncells
    nn = 1
    do n = 1, maxval(tlength)
      tdistr(n) = count(tlength == n)
      cummt = cummt + tdistr(n)
      if (cummt >= r_bucketsize(nn)) then
        bucketsize(nn) = cummt
        bucket_max(nn) = n
        nn = nn + 1
      end if
    end do
    bucket_min(1)          = 1
    bucket_min(2:nbuckets) = bucket_max(1:nbuckets-1) + 1

    do n = 2, nbuckets
      bucketsize(n) = bucketsize(n) - sum(bucketsize(1:n-1))
    end do
    bucketsize = max(bucketsize,0)
!Accumulate the individual cell statistics
    do n = 1, nbuckets
      write (*,*) 'Process bucket', n
      if (bucketsize(n) > 0) then
        allocate(base(bucket_max(n), bucketsize(n)))
        base = huge(1.)
        allocate(top (bucket_max(n), bucketsize(n)))
        top  = fillvalue_r
        allocate(area(bucket_max(n), bucketsize(n)))
        area = 0.
        allocate(maxarea(bucket_max(n), bucketsize(n)))
        maxarea = fillvalue_r
        allocate(maxarealoc(bucket_max(n), bucketsize(n)))
        maxarealoc = fillvalue_r
        allocate(recon(ceiling(maxheight/dz), bucket_max(n)))
        allocate(vol (bucket_max(n), bucketsize(n)))
        vol  = 0.
        allocate(val (bucket_max(n), bucketsize(n)))
        val  = 0.
        allocate(xcenter (bucket_max(n), bucketsize(n)))
        xcenter = fillvalue_r
        allocate(ycenter (bucket_max(n), bucketsize(n)))
        ycenter = fillvalue_r
        allocate(icenter(bucket_max(n)))
        allocate(jcenter(bucket_max(n)))
        allocate(ianchor(bucket_max(n)))
        allocate(janchor(bucket_max(n)))
        allocate(npts   (bucket_max(n)))

        allocate(duration(bucketsize(n)))
        allocate(mintime(bucketsize(n)))
        allocate(maxtime(bucketsize(n)))
        allocate(id(bucketsize(n)))

        nn   = 0
        iret = firstcell(cell)
        do
          if (iret == -1) exit
          tmin = minval(cell%loc(3,1:cell%nelements))
          tmax = maxval(cell%loc(3,1:cell%nelements))
          if (tmax - tmin + 1 >= bucket_min(n) .and. tmax - tmin + 1 <= bucket_max(n)) then
            nn = nn + 1
            if (mod(nn,1000)==0) then
              write (*,'(A,i10,A,i10)') '..Cell nr', nn,'/',bucketsize(n)
            end if
            id      (nn)  = cell%id
            duration(nn) = dt * (tmax - tmin)
            mintime(nn)  = dt * tmin
            maxtime(nn)  = dt * tmax

            base    (tmax-tmin+2:bucket_max(n), nn) = fillvalue_r
            top     (tmax-tmin+2:bucket_max(n), nn) = fillvalue_r
            area    (tmax-tmin+2:bucket_max(n), nn) = fillvalue_r
            vol     (tmax-tmin+2:bucket_max(n), nn) = fillvalue_r
            val     (tmax-tmin+2:bucket_max(n), nn) = fillvalue_r
            recon = 0.
!             icenter = 0
!             jcenter = 0
!             npts  = 0

            do nel = 1, cell%nelements
              i = cell%loc(1,nel)
              j = cell%loc(2,nel)
              t = cell%loc(3,nel)
              tt = t - tmin + 1
              if (cell%id<0) then
                write (*,*) 'DANGER: Cell id < 0'
              else
                slab(i,j,t) = cell%id
              end if
              rbase = real(cell%value(ibase,nel)) * heightrange + heightzero
              rtop  = real(cell%value(itop,nel)) * heightrange + heightzero
              base(tt, nn) = min(base(tt, nn), rbase)
              top (tt, nn) = max(top (tt, nn), rtop)
              vol (tt, nn) = vol(tt, nn) + rtop-rbase
              area(tt, nn) = area(tt, nn) + 1.
              recon(floor(rbase/dz)+1:floor(rtop/dz)+1,tt) = recon(floor(rbase/dz)+1:floor(rtop/dz)+1,tt) + 1.
              val (tt, nn) = val(tt, nn) +  real(cell%value(ivalue,nel))
              !Calculate the center of the cloud
              if (npts(tt) == 0) then
                ianchor(tt) = i
                janchor(tt) = j
              end if
              npts(tt) = npts(tt) + 1
              idev = i-ianchor(tt)
              icenter(tt) = icenter(tt) + idev
              if (abs(idev) > nx/2) then
                icenter(tt) = icenter(tt) - sign(nx,idev)
              end if
              jdev = j-janchor(tt)
              jcenter(tt) = jcenter(tt) + jdev
              if (abs(jdev) > ny/2) then
                jcenter(tt) = jcenter(tt) - sign(ny,jdev)
              end if
            end do
            do tt = 1, tmax-tmin+1
              xcenter(tt,nn) = (real(ianchor(tt))+real(icenter(tt))/real(npts(tt)) - 0.5*(nx-1))*dx
              if (xcenter(tt,nn) < -0.5*real(nx)*dx) then
                xcenter(tt,nn) = xcenter(tt,nn)+real(nx)*dx
              elseif (xcenter(tt,nn) > 0.5*real(nx)*dx) then
                xcenter(tt,nn) = xcenter(tt,nn) - real(nx)*dx
              end if

              ycenter(tt,nn) = (real(janchor(tt))+real(jcenter(tt))/real(npts(tt)) - 0.5*(ny-1))*dy
              if (ycenter(tt,nn) < -0.5*real(ny)*dy) then
                ycenter(tt,nn) = ycenter(tt,nn)+real(ny)*dy
              elseif (ycenter(tt,nn) > 0.5*real(ny)*dy) then
                ycenter(tt,nn) = ycenter(tt,nn) - real(ny)*dy
              end if
              maxarea(tt,nn) = maxval(recon(:,tt)) *dx*dy
              maxarealoc(tt,nn) = maxloc(recon(:,tt),1)*dz
            end do
          end if
          iret = nextcell(cell)
        end do
        where(area>0)
          val = val/area*valrange+valzero
          area = area * dx * dy
          vol  = vol * dx * dy
        elsewhere
          val = fillvalue_r
          area = fillvalue_r
          vol = fillvalue_r
          base = fillvalue_r
        end where
        write (*,*) '..write to disk'
        !Write output to netcdf
        !Write nr dimension for the current block
        nrnc%name     = trim(bucketname(n))//trim(ovarstem%name)
        nrnc%longname = trim(bucketlongname(n))//' '//trim(ovarstem%longname)//'s'
        nrnc%units    = '#'
        nrnc%dim(1) = bucketsize(n)
        call define_ncdim(fid, nrnc, nf90_int)
        call write_ncvar(fid, nrnc, (/(i, i=1,nrnc%dim(1))/))

        !Write time dimension for the current block
        timenc%name     = trim(bucketname(n))//trim(ovarstem%name)//'t'
        timenc%longname = trim(bucketlongname(n))//' '//trim(ovarstem%name)// ' time'
        timenc%units    = tnc%units
        timenc%dim(1) = bucket_max(n)
        call define_ncdim(fid, timenc)
        call write_ncvar(fid, timenc, dt*(/(i, i=1,timenc%dim(1))/))

        !Set the general netcdf settings of the variables only depending on cellnr
        allocate(ovar%dim(1))
        ovar%dim      = shape(id)
        allocate(ovar%dimids(1))
        ovar%dimids(1) = nrnc%dimids(1)

        !Write to netcdf file: Cell id
        ovar%name     = trim(nrnc%name)//'id'
        ovar%longname = trim(nrnc%longname)//' id'
        ovar%units    = '#'
        call define_ncvar(fid, ovar, nf90_int)
        call write_ncvar(fid, ovar, id)

        !Write to netcdf file: Cell duration
        ovar%name     = trim(nrnc%name)//'dur'
        ovar%longname = trim(nrnc%longname)//' duration'
        ovar%units    = tnc%units
        call define_ncvar(fid, ovar, nf90_float)
        call write_ncvar(fid, ovar, duration)

        !Write to netcdf file: Cell mintime
        ovar%name     = trim(nrnc%name)//'tmin'
        ovar%longname = trim(nrnc%longname)//' time of appearance'
        ovar%units    = tnc%units
        call define_ncvar(fid, ovar, nf90_float)
        call write_ncvar(fid, ovar, mintime)

        !Write to netcdf file: Cell maxtime
        ovar%name     = trim(nrnc%name)//'tmax'
        ovar%longname = trim(nrnc%longname)//' time of extinction'
        ovar%units    = tnc%units
        call define_ncvar(fid, ovar, nf90_float)
        call write_ncvar(fid, ovar, maxtime)

        deallocate(ovar%dim, ovar%dimids)
        !Set the general netcdf settings of the variables only depending on cellnr and time
        allocate(ovar%dim(2))
        ovar%dim      = shape(base)
        allocate(ovar%dimids(2))
        ovar%dimids(1) = timenc%dimids(1)
        ovar%dimids(2) = nrnc%dimids(1)

        !Write to netcdf file: Cell base
        ovar%name     = trim(nrnc%name)//'x'
        ovar%longname = trim(nrnc%longname)//' x coordinate of the center of mass'
        ovar%units    = 'm'
        call define_ncvar(fid, ovar, nf90_float)
        call write_ncvar(fid, ovar, xcenter)

        !Write to netcdf file: Cell base
        ovar%name     = trim(nrnc%name)//'y'
        ovar%longname = trim(nrnc%longname)//' y coordinate of the center of mass'
        ovar%units    = 'm'
        call define_ncvar(fid, ovar, nf90_float)
        call write_ncvar(fid, ovar, ycenter)

        !Write to netcdf file: Cell base
        ovar%name     = trim(nrnc%name)//'base'
        ovar%longname = trim(nrnc%longname)//' base'
        ovar%units    = 'm'
        call define_ncvar(fid, ovar, nf90_float)
        call write_ncvar(fid, ovar, base)

        !Write to netcdf file: Cell top
        ovar%name     = trim(nrnc%name)//'top'
        ovar%longname = trim(nrnc%longname)//' top'
        ovar%units    = 'm'
        call define_ncvar(fid, ovar, nf90_float)
        call write_ncvar(fid, ovar, top)

        !Write to netcdf file: Cell area
        ovar%name     = trim(nrnc%name)//'area'
        ovar%longname = trim(nrnc%longname)//' Area'
        ovar%units    = 'm^2'
        call define_ncvar(fid, ovar, nf90_float)
        call write_ncvar(fid, ovar, area)

        !Write to netcdf file: Cell maxarea
        ovar%name     = trim(nrnc%name)//'maxarea'
        ovar%longname = trim(nrnc%longname)//' Maximum Area'
        ovar%units    = 'm^2'
        call define_ncvar(fid, ovar, nf90_float)
        call write_ncvar(fid, ovar, maxarea)

        !Write to netcdf file: Cell max area loc
        ovar%name     = trim(nrnc%name)//'maxarealoc'
        ovar%longname = trim(nrnc%longname)//' Location of the maximum area'
        ovar%units    = 'm'
        call define_ncvar(fid, ovar, nf90_float)
        call write_ncvar(fid, ovar, maxarealoc)

        !Write to netcdf file: Cell vol
        ovar%name     = trim(nrnc%name)//'vol'
        ovar%longname = trim(nrnc%longname)//' Volume'
        ovar%units    = 'm^3'
        call define_ncvar(fid, ovar, nf90_float)
        call write_ncvar(fid, ovar, vol)

        !Write to netcdf file: Cell val
        ovar%name     = trim(nrnc%name)//trim(ivar(icell)%name)
        ovar%longname = trim(nrnc%longname)//' mean '//trim(ivar(icell)%longname)
        ovar%units    = ivar(icell)%units
        call define_ncvar(fid, ovar, nf90_float)
        call write_ncvar(fid, ovar, val)

        deallocate(ovar%dim, ovar%dimids)

        deallocate(base, top, area, maxarea, maxarealoc, recon, vol, val)
        deallocate(icenter,jcenter, xcenter, ycenter,ianchor, janchor, npts)
        deallocate(duration, mintime, maxtime, id)



        write (*,*) '..Relationships'
        allocate(nrelatives(bucketsize(n)))
        allocate(relatives(nrel_max,bucketsize(n)))
        relatives = fillvalue_i


        nn   = 0
        iret = firstcell(cell)
        do
          if (iret == -1) exit
          tmin = minval(cell%loc(3,1:cell%nelements))
          tmax = maxval(cell%loc(3,1:cell%nelements))
          if (tmax - tmin + 1 >= bucket_min(n) .and. tmax - tmin + 1 <= bucket_max(n)) then
            nn = nn + 1
!             if (mod(nn,1000)==0) then
!               write (*,'(A,i10,A,i10)') '..Cell nr', nn,'/',bucketsize(n)
!             end if

            nrelatives(nn)  = cell%nparents
            do nel = 1, cell%nparents
              relatives(nel, nn) = cell%parents(nel)%p%id
            end do
          end if
          iret = nextcell(cell)
        end do
        !Write to netcdf file: Nr parents
        if (any(nrelatives>0)) then
          allocate(ovar%dim(1))
          ovar%dim      = shape(nrelatives)
          allocate(ovar%dimids(1))
          ovar%dimids(1) = nrnc%dimids(1)
          ovar%name     = trim(nrnc%name)//'nrpar'
          ovar%longname = trim(nrnc%longname)//' number of parents'
          ovar%units    = '#'
          call define_ncvar(fid, ovar, nf90_int)
          call write_ncvar(fid, ovar, nrelatives)
          deallocate(ovar%dim, ovar%dimids)
        end if
        !Write to netcdf file: Parents
        if (any(nrelatives>0)) then
          allocate(ovar%dim(2))
          ovar%dim      = shape(relatives)
          allocate(ovar%dimids(2))
          ovar%dimids(1) = relnc%dimids(1)
          ovar%dimids(2) = nrnc%dimids(1)
          ovar%name     = trim(nrnc%name)//'par'
          ovar%longname = trim(nrnc%longname)//' parent ids'
          ovar%units    = '#'
          call define_ncvar(fid, ovar, nf90_int)
          call write_ncvar(fid, ovar, relatives)
          deallocate(ovar%dim, ovar%dimids)
        end if

        !Write to netcdf file: Nr children
        relatives = fillvalue_i
        nn   = 0
        iret = firstcell(cell)
        do
          if (iret == -1) exit
          tmin = minval(cell%loc(3,1:cell%nelements))
          tmax = maxval(cell%loc(3,1:cell%nelements))
          if (tmax - tmin + 1 >= bucket_min(n) .and. tmax - tmin + 1 <= bucket_max(n)) then
            nn = nn + 1
            nrelatives(nn)  = cell%nchildren
            do nel = 1, cell%nchildren
              relatives(nel, nn) = cell%children(nel)%p%id
            end do
          end if
          iret = nextcell(cell)
        end do
        if (any(nrelatives>0)) then
          allocate(ovar%dim(1))
          ovar%dim      = shape(nrelatives)
          allocate(ovar%dimids(1))
          ovar%dimids(1) = nrnc%dimids(1)
          ovar%name     = trim(nrnc%name)//'nrchild'
          ovar%longname = trim(nrnc%longname)//' number of children'
          ovar%units    = '#'
          call define_ncvar(fid, ovar, nf90_int)
          call write_ncvar(fid, ovar, nrelatives)
          deallocate(ovar%dim, ovar%dimids)
        end if
        !Write to netcdf file: Children
        if (any(nrelatives>0)) then
          allocate(ovar%dim(2))
          ovar%dim      = shape(relatives)
          allocate(ovar%dimids(2))
          ovar%dimids(1) = relnc%dimids(1)
          ovar%dimids(2) = nrnc%dimids(1)
          ovar%name     = trim(nrnc%name)//'child'
          ovar%longname = trim(nrnc%longname)//' child ids'
          ovar%units    = '#'
          call define_ncvar(fid, ovar, nf90_int)
          call write_ncvar(fid, ovar, relatives)
          deallocate(ovar%dim, ovar%dimids)
        end if
!
!         !Write to netcdf file: Nr siblings
! !         if (lsiblings) then
          nrelatives = fillvalue_i
          nn   = 0
          iret = firstcell(cell)
!
          do
            if (iret == -1) exit
            tmin = minval(cell%loc(3,1:cell%nelements))
            tmax = maxval(cell%loc(3,1:cell%nelements))
            if (tmax - tmin + 1 >= bucket_min(n) .and. tmax - tmin + 1 <= bucket_max(n)) then
              nn = nn + 1
              nrelatives(nn)  = cell%cloudtype
! !               do nel = 1, cell%nsiblings
! !                 relatives(nel, nn) = cell%siblings(nel)%p%id
! !               end do
            end if
            iret = nextcell(cell)
          end do
          if (any(nrelatives>0)) then
            allocate(ovar%dim(1))
            ovar%dim      = shape(nrelatives)
            allocate(ovar%dimids(1))
            ovar%dimids(1) = nrnc%dimids(1)
            ovar%name     = trim(nrnc%name)//'type'
            ovar%longname = trim(nrnc%longname)//' type. 1 for passive, 2 for single pulse, 3 for outflow, 4 for active'
            ovar%units    = '#'
            call define_ncvar(fid, ovar, nf90_int)
            call write_ncvar(fid, ovar, nrelatives)
            deallocate(ovar%dim, ovar%dimids)
          end if

!         !Write to netcdf file: Cloudsystemid
          nrelatives = fillvalue_i
          nn   = 0
          iret = firstcell(cell)
          do
            if (iret == -1) exit
            tmin = minval(cell%loc(3,1:cell%nelements))
            tmax = maxval(cell%loc(3,1:cell%nelements))
            if (tmax - tmin + 1 >= bucket_min(n) .and. tmax - tmin + 1 <= bucket_max(n)) then
              nn = nn + 1
              nrelatives(nn)  = cell%cloudsystemnr
            end if
            iret = nextcell(cell)
          end do
          if (any(nrelatives>0)) then
            allocate(ovar%dim(1))
            ovar%dim      = shape(nrelatives)
            allocate(ovar%dimids(1))
            ovar%dimids(1) = nrnc%dimids(1)
            ovar%name     = trim(nrnc%name)//'sysid'
            ovar%longname = trim(nrnc%longname)//' id of the system'
            ovar%units    = '#'
            call define_ncvar(fid, ovar, nf90_int)
            call write_ncvar(fid, ovar, nrelatives)
            deallocate(ovar%dim, ovar%dimids)
          end if
!           !Write to netcdf file: siblings
!           if (any(nrelatives>0)) then
!             allocate(ovar%dim(2))
!             ovar%dim      = shape(relatives)
!             allocate(ovar%dimids(2))
!             ovar%dimids(1) = relnc%dimids(1)
!             ovar%dimids(2) = nrnc%dimids(1)
!             ovar%name     = trim(nrnc%name)//'sib'
!             ovar%longname = trim(nrnc%longname)//' sibling ids'
!             ovar%units    = '#'
!             call define_ncvar(fid, ovar, nf90_int)
!             call write_ncvar(fid, ovar, relatives)
!             deallocate(ovar%dim, ovar%dimids)
!           end if
!         end if
        !Write to netcdf file: Nr splitters
        nrelatives = 0
        relatives = fillvalue_i
        nn   = 0
        iret = firstcell(cell)
        do
          if (iret == -1) exit
          tmin = minval(cell%loc(3,1:cell%nelements))
          tmax = maxval(cell%loc(3,1:cell%nelements))
          if (tmax - tmin + 1 >= bucket_min(n) .and. tmax - tmin + 1 <= bucket_max(n)) then
            nn = nn + 1
            nrelatives(nn)  = cell%nsplitters
            do nel = 1, cell%nsplitters
              relatives(nel, nn) = cell%splitters(nel)%p%id
            end do
          end if
          iret = nextcell(cell)
        end do
        if (any(nrelatives>0)) then
          allocate(ovar%dim(1))
          ovar%dim      = shape(nrelatives)
          allocate(ovar%dimids(1))
          ovar%dimids(1) = nrnc%dimids(1)
          ovar%name     = trim(nrnc%name)//'nrsplit'
          ovar%longname = trim(nrnc%longname)//' number of splitters'
          ovar%units    = '#'
          call define_ncvar(fid, ovar, nf90_int)
          call write_ncvar(fid, ovar, nrelatives)
          deallocate(ovar%dim, ovar%dimids)
        end if
        !Write to netcdf file: splitters
        if (any(nrelatives>0)) then
          allocate(ovar%dim(2))
          ovar%dim      = shape(relatives)
          allocate(ovar%dimids(2))
          ovar%dimids(1) = relnc%dimids(1)
          ovar%dimids(2) = nrnc%dimids(1)
          ovar%name     = trim(nrnc%name)//'split'
          ovar%longname = trim(nrnc%longname)//' splitter ids'
          ovar%units    = '#'
          call define_ncvar(fid, ovar, nf90_int)
          call write_ncvar(fid, ovar, relatives)
          deallocate(ovar%dim, ovar%dimids)
        end if
        deallocate(nrelatives, relatives)


      end if
    end do

    !write slab
    write (*,*) '..Write slab'
    allocate(ovar%dim(3))
    ovar%dim      = shape(slab)
    allocate(ovar%dimids(3))
    ovar%dimids(1) = xnc%dimids(1)
    ovar%dimids(2) = ync%dimids(1)
    ovar%dimids(3) = tnc%dimids(1)
    ovar%name = 'nr'//trim(ovarstem%name)
    ovar%longname = trim(ovarstem%name) // ' number'
    ovar%units = '#'
    call define_ncvar(fid, ovar, nf90_int)
    call write_ncvar(fid, ovar, slab)
    write (*,*) '.. leaving statistics'
  end subroutine dostatistics

end module modtrack

!> Main entrypoint of cloud tracking code
!>
program tracking
  use modnetcdf
  use modtrack

  implicit none

  real, parameter :: thermmin = -1.
  real, parameter :: thermmax = 10000.
  real, parameter :: thermthres = 300.

  real, parameter :: lwpmin   = -1.
  real, parameter :: lwpmax   = 10.
  real, parameter :: lwpthres = 0.01

  real, parameter :: coremin  = -5.
  real, parameter :: coremax  = 5.
  real, parameter :: corethres = 0.5

  real, parameter :: rwpmin   = -1.
  real, parameter :: rwpmax   = 10.
  real, parameter :: rwpthres = 0.01

  real, parameter :: distmin = -1.
  real, parameter :: distmax = 5000.

  real, parameter :: maxheight = 5000.

  integer, parameter :: nchunk = 100

  integer, parameter :: nmincells_cloud  = 1
  integer, parameter :: nmincells        = 4

  ! center and range values for parameters
  real, parameter :: thermzero  = 0.5*(thermmax + thermmin)
  real, parameter :: thermrange = (thermmax - thermmin)/real(huge(1_2))

  real, parameter :: lwpzero  = 0.5*(lwpmax + lwpmin)
  real, parameter :: lwprange = (lwpmax - lwpmin)/real(huge(1_2))

  real, parameter :: corezero  = 0.5*(coremax + coremin)
  real, parameter :: corerange = (coremax - coremin)/real(huge(1_2))

  real, parameter :: rwpzero  = 0.5*(rwpmax + rwpmin)
  real, parameter :: rwprange = (rwpmax - rwpmin)/real(huge(1_2))

  real, parameter :: distzero  = 0.5*(distmax + distmin)
  real, parameter :: distrange = (distmax - distmin)/real(huge(1_2))

  integer(kind=2), parameter :: i_lwpthres   = (lwpthres - lwpzero)/lwprange
  integer(kind=2), parameter :: i_corethres  = (corethres - corezero)/corerange
  integer(kind=2), parameter :: i_thermthres = (thermthres - thermzero)/thermrange
  integer(kind=2), parameter :: i_rwpthres   = (rwpthres - rwpzero)/rwprange

  ! runtime variables
  type(cellptr), dimension(:,:,:), allocatable :: parentarr

  logical :: lcore = .false., lcloud = .false., lthermal = .false., lrain = .false., lrwp = .false.
  integer :: i,j,k,n, kk, kkmax
  integer :: ub, lb, vlength, fid, finput, finput2
  character(100) :: criterion, filename, stem, ctmp
  real    :: heightmin, heightrange, valmin, valrange
  type(netcdfvar), dimension(10) :: ivar
  type(netcdfvar) :: ovar
  real, allocatable, dimension(:) :: x, y, t
  integer(kind=2), dimension(:,:,:), allocatable :: base, top
  integer(kind=2), dimension(:), allocatable :: minbasecloud, minbasetherm
  real, dimension(:,:,:), allocatable :: readfield, readfield2


  if (command_argument_count() == 0) then
     print *, "./tracking [filename-base] [starting timestep] [final timestep] [analysis variables]"
     print *, ""
     print *, "analysis variables: `core`, `cloud`, `liquid`??, `thermal`, `rain` and `all`"
     print *, "analysis variables defines which fields will be analysed"
     call exit(-1)
  endif

  ! defined in `modtrack`
  nt = 0
  tstart = 1
  cbstep = (300.)/distrange
!  cbstep = (3000.)/distrange  !AXEL

  call parse_command_arguments()
  call read_input(stem)

  allocate(var(nx, ny, tstart:nt, nvar))
  allocate(parentarr(nx, ny, tstart:nt))
  allocate(bool(nx, ny, tstart:nt))
  allocate(readfield(nx, ny, nchunk))
  if (lrwp) allocate(readfield2(nx, ny, nchunk))

  call setup_output(stem)


  minparentel = 100!nint(50./(dx*dy*dt))
  if (lthermal) then
    write (*,*) 'Thermals....'
    ivar(1)%name  = 'trcbase'
    ivar(2)%name  = 'trctop'

    call read_named_input('trcpath', -thermzero, thermrange, thermmin, fillvalue_i16, ithermal)

    call read_named_input('trcbase', -distzero, distrange, distmin, fillvalue_i16, ibase)
    call read_named_input('trctop', -distzero, distrange, distmin, fillvalue_i16, itop)

    allocate(minbasetherm(tstart:nt))
    minbasetherm = (100.-distzero)/distrange

    !Loop over thermals
    bool = -1
    where (var(:,:,:,ivalue) > i_thermthres)
      bool = 0
    end where
!     do n=1,tstart-1
!       bool(:,:,n) = -1
!     end do
!       call checkframes
    call dotracking(thermal, nthermals,nmincells)
  end if

  if (lcore .or. lcloud) then
    ivar(1)%name  = 'cldbase'
    ivar(2)%name  = 'cldtop'

    call read_named_input('lwp', -lwpzero, lwprange, lwpmin, fillvalue_i16, ivalue)

    call read_named_input('cldbase', -distzero, distrange, distmin, fillvalue_i16, ibase)
    call read_named_input('cldtop', -distzero, distrange, distmin, fillvalue_i16, itop)


    !Loop over cloud cores
    if (lcore) then
    ! Read in core specific input
      write (*,*) 'Cores....'
      ivalue = 4
      call read_named_input('core', -corezero, corerange, coremin, fillvalue_i16, ivalue)

      bool = -1
      allocate(minbasecloud(tstart:nt))

      do n=tstart, nt
        minbasecloud(n) = (maxval(var(:,:,n,itop)) + minval(var(:,:,n,ibase),var(:,:,n,ibase)>fillvalue_i16) )/2
        if (any(var(:,:,n,3) > i_lwpthres)) then
          where (var(:,:,n,3) > i_lwpthres &
         .and. var(:,:,n,ivalue) > i_corethres &
         .and. var(:,:,n,ibase) < minbasecloud(n))
            bool(:,:,n) = 0
          end where
        end if
      end do

!       call checkframes
      call dotracking(core, ncores,nmincells)
      ivalue = 3
   end if
    !Loop over clouds
    if (lcloud) then
    write (*,*) 'Clouds....'
      bool = -1
      where (var(:,:,:,ilwp) > i_lwpthres)
        bool = 0
      end where
!       do n=1,tstart-1
!         bool(:,:,n) = -1
!       end do
!       call checkframes
      if (lcore) then
        call fillparentarr(core, minbasecloud, parentarr)
        call dotracking(cloud, nclouds,nmincells_cloud, parentarr)
      else
        call dotracking(cloud, nclouds,nmincells_cloud)
      end if
      if (lthermal) then
        allocate(base(nx, ny, tstart:nt))
        allocate(top(nx, ny, tstart:nt))
        call fillparentarr(thermal, minbasetherm, parentarr, base, top)
        call findparents(cloud, parentarr, base, top)
        deallocate(base,top)
      end if
    end if
  end if

  if (lrain) then
    write (*,*) 'Rain....'
    filename = trim(stem)//trim(ivar(irain)%name)//'.nc'
    write(*,*) 'Reading ', trim(filename)
    call check ( nf90_open (trim(filename), NF90_NOWRITE, finput) )
    call inquire_ncvar(finput, ivar(irain))

    do k = tstart,nt,nchunk
      write (*,*) 'Reading t = ',k
      kkmax = min(nt-k+1,nchunk)
      call read_ncvar(finput, ivar(irain), readfield,(/1,1,k/),(/nx,ny,kkmax/))
      do kk = 1,kkmax
        do j = 1, ny
          do i = 1, nx
            if (readfield(i,j,kk) >= rwpmin) then
              var(i,j,k+kk-1,ivalue) = (readfield(i,j,kk) - rwpzero ) / rwprange
            else
              var(i,j,k+kk-1,ivalue) = fillvalue_i16
            end if
          end do
        end do
      end do
    end do
!     call check(nf90_close(finput))
    ivar(1)%name  = 'rwpbase'
    ivar(2)%name  = 'rwptop'

    do n = 1,2
      filename = trim(stem)//trim(ivar(n)%name)//'.nc'
      write(*,*) 'Reading ', trim(filename)
      call check ( nf90_open (trim(filename), NF90_NOWRITE, finput) )
      call inquire_ncvar(finput, ivar(n))
      do k = tstart,nt,nchunk
        write (*,*) 'Reading t = ',k
        kkmax = min(nt-k+1,nchunk)
!         call read_ncvar(finput, ivar(n), readfield,(/1,1,k/),(/nx,ny,kkmax/))
        do kk = 1,kkmax
          do j = 1, ny
            do i = 1, nx
!               if (readfield(i,j,kk) >= distmin) then
!                 var(i,j,k+kk-1,ivalue) = (readfield(i,j,kk) - distzero ) / distrange
!               else
!                 var(i,j,k+kk-1,ivalue) = fillvalue_i16
!               end if
            var(i,j,k+kk-1,1) = (0. - distzero ) / distrange
            var(i,j,k+kk-1,2) = (4000. - distzero ) / distrange
             end do
          end do
        end do
      end do
    end do

!     deallocate(readfield)

    !Loop over rain patches
    bool = -1
    where (var(:,:,:,ivalue) > i_rwpthres)
      bool = 0
    end where
!       do n=1,tstart-1
!         bool(:,:,n) = -1
!       end do
!       call checkframes
    call dotracking(rain, nrains,nmincells)
    if (lcloud) then
      allocate(base(nx, ny, tstart:nt))
      allocate(top(nx, ny, tstart:nt))
      if (.not. allocated(minbasecloud)) allocate(minbasecloud(tstart:nt))
      minbasecloud = huge(1_2)
      call fillparentarr(cloud, minbasecloud, parentarr, base, top)
      call findparents(rain, parentarr, base, top)
      deallocate(base,top)
    end if
  end if

  deallocate(readfield, parentarr, var, bool)

  if (lthermal) then
    ovar%name     = 'thrm'
    ovar%longname = 'Thermal'
    call dostatistics(thermal, nthermals, ithermal, fid, ivar, distzero, distrange, maxheight,thermzero, thermrange, ovar)
!     call delete_all(thermal)

  end if
  if (lcore) then
    ovar%name     = 'core'
    ovar%longname = 'Core'
    call dostatistics(core, ncores, icore, fid, ivar, distzero, distrange, maxheight,corezero, corerange, ovar)
!     call delete_all(core)
  end if
  if (lcloud) then
    ovar%name     = 'cloud'
    ovar%longname = 'Cloud'
    call dostatistics(cloud, nclouds, ilwp, fid, ivar, distzero, distrange, maxheight,lwpzero, lwprange, ovar)
!     call delete_all(cloud)
  end if
  if (lrain) then
    ovar%name     = 'rain'
    ovar%longname = 'Rain patch'
    call dostatistics(rain, nrains, irain, fid, ivar, distzero, distrange, maxheight, rwpzero, rwprange, ovar)
!     call delete_all(rain)
  end if
  call check ( nf90_close(fid))
  write (*,*) '..Done'

  contains
     subroutine setup_output(stem)
        character(len=100), intent(in) :: stem

        call check ( nf90_create(trim(stem)//'track.nc', NF90_HDF5, fid))
        if (lthermal) then
           call check ( nf90_put_att(fid, nf90_global, 'Thermal threshold',thermthres))
           call check ( nf90_put_att(fid, nf90_global, 'Thermal Min size',nmincells))
        end if
        if (lcloud) then
           call check ( nf90_put_att(fid, nf90_global, 'LWP threshold',lwpthres))
           call check ( nf90_put_att(fid, nf90_global, 'LWP Min size',nmincells_cloud))
        end if
        if (lcore) then
           call check ( nf90_put_att(fid, nf90_global, 'Core threshold',corethres))
           call check ( nf90_put_att(fid, nf90_global, 'Core Min size',nmincells))
        end if
        if (lrain) then
           call check ( nf90_put_att(fid, nf90_global, 'RWP threshold',rwpthres))
           call check ( nf90_put_att(fid, nf90_global, 'RWP Min size',nmincells))
        end if
        call check ( nf90_enddef(fid))
     end subroutine setup_output

     subroutine read_input(stem)
        character(len=100), intent(in) :: stem

        filename = trim(stem)//trim(ivar(nvar)%name)//'.nc'
        write(*,*) 'Reading dimensions, filename=', filename
        call check ( nf90_open (trim(filename), NF90_NOWRITE, finput) )
        tnc%name = 'time'
        call inquire_ncvar(finput, tnc)
        if (nt > 0) then
           tnc%dim(1) = nt -tstart + 1
        else
           nt = tnc%dim(1)
        end if
        allocate(t(tstart:nt))
        call read_ncvar(finput, tnc, t,(/tstart/),(/nt-tstart+1/))
        call define_ncdim(fid, tnc)
        call write_ncvar(fid, tnc, t)
        dt = t(tstart+1)-t(tstart)

        xnc%name = 'xt'
        call inquire_ncvar(finput, xnc)
        nx = xnc%dim(1)
        allocate(x(nx))
        call read_ncvar(finput, xnc, x)
        call define_ncdim(fid, xnc)
        call write_ncvar(fid, xnc, x)
        dx = x(2)-x(1)

        ync%name = 'yt'
        call inquire_ncvar(finput, ync)
        ny = ync%dim(1)
        allocate(y(ny))
        call read_ncvar(finput, ync, y)
        call define_ncdim(fid, ync)
        call write_ncvar(fid, ync, y)
        dy = y(2)-y(1)
        if (lcore) then
           nvar = 4
        else
           nvar = 3
        end if
        call check(nf90_close(finput))
        write (*,*) nx, ny, nt, tstart,nvar
     end subroutine read_input

     subroutine read_named_input(var_name, value_offset, value_scaling, min_value, fill_value, var_index)
        character(len=*), intent(in) :: var_name
        integer, intent(in) :: var_index
        real, intent(in) :: value_offset, value_scaling, min_value
        integer(kind=2) :: fill_value

        character(len=200) :: filename
        type(netcdfvar) :: nc_var
        real, dimension(:,:,:), allocatable :: var_data

        nc_var%name = var_name
        filename = trim(stem)//trim(var_name)//'.nc'

        write(*,*) 'Reading ', trim(filename)
        call check ( nf90_open (trim(filename), NF90_NOWRITE, finput) )
        call inquire_ncvar(finput, nc_var)
        do k = tstart,nt,nchunk
          write (*,*) 'Reading t = ',k
          kkmax = min(nt-k+1,nchunk)
          call read_ncvar(finput, nc_var, var_data,(/1,1,k/),(/nx,ny,kkmax/))
          do kk = 1,kkmax
            do j = 1, ny
              do i = 1, nx
                if (var_data(i,j,kk) >= min_value) then
                  var(i,j,k+kk-1,var_index) = (var_data(i,j,kk) + value_offset ) / value_scaling
                else
                  var(i,j,k+kk-1,var_index) = fill_value
                end if
              end do
            end do
          end do
        end do
        call check(nf90_close(finput))
     end subroutine read_named_input

     ! Parse command line arguments
     ! `./tracking_time [filename-base] [starting timestep] [final timestep] [analysis variables]
     !
     ! analysis variables: `core`, `cloud`, `liquid`??, `thermal`, `rain` and `all`
     ! analysis variables defines which fields will be analysed
     subroutine parse_command_arguments()
        call get_command_argument(1,stem)
        stem = trim(stem)//'.out.xy.'
        call get_command_argument(2,ctmp)
        read(ctmp,'(i4)') tstart
        call get_command_argument(3,ctmp)
        read(ctmp,'(i4)') nt
        call get_command_argument(4,criterion)

        nvar = 2
        nrel_max = 0
        vlength = len_trim(criterion)
        lb = 1
        ub = -1
        n=0
        do
        if (ub == vlength) exit
        lb = ub + 2
        ub = index(criterion(lb:vlength),',') - 2 + lb !Comma separated list of variables
        if (ub < lb) ub = vlength
        select case (trim(adjustl(criterion(lb:ub))))
        case ('core')
           lcore = .true.
           if (ilwp < 0) then
              nvar = nvar + 1
              ilwp = nvar
           end if
           if (icore < 0) then
              nvar = nvar + 1
              icore = nvar
           end if

        case ('cloud')
           lcloud = .true.
           if (ilwp < 0) then
              nvar = nvar + 1
              ilwp = nvar
           end if
        case ('liquid')
           lcloud = .true.
           lrwp   = .true.
           if (ilwp < 0) then
              nvar = nvar + 1
              ilwp = nvar
           end if
           if (irain < 0) then
              nvar  = nvar + 1
              irain = nvar
           end if
        case ('thermal')
           lthermal = .true.
           if (ithermal < 0) then
              nvar = nvar + 1
              ithermal = nvar
           end if
        case ('rain')
           lrain = .true.
           if (irain < 0) then
              nvar = nvar + 1
              irain = nvar
           end if
        case ('all')
           lcore = .true.
           lcloud = .true.
           lthermal = .true.
           lrain = .true.
           if (ilwp < 0) then
              nvar = nvar + 1
              ilwp = nvar
           end if
           if (icore < 0) then
              nvar = nvar + 1
              icore = nvar
           end if
           if (ithermal < 0) then
              nvar = nvar + 1
              ithermal = nvar
           end if
           if (irain < 0) then
              nvar = nvar + 1
              irain = nvar
           end if
        end select
        end do
        if (ilwp     > 0) ivar(ilwp)%name     = 'lwp'
        if (icore    > 0) ivar(icore)%name    = 'core'
        if (ithermal > 0) ivar(ithermal)%name = 'trcpath'
        if (irain    > 0) ivar(irain)%name    = 'rwp'
        ! End parse command line arguments
     end subroutine parse_command_arguments
end program tracking
