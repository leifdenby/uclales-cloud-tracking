module modtrack
  use tracking_common, only: celltype, cellptr

  use tracking_common, only: dt, dx, dy
  use tracking_common, only: ibase, itop, ivalue
  use tracking_common, only: nrel_max, nt, nx, ny
  use tracking_common, only: tstart
  use tracking_common, only: createcell, deletecell, firstcell, nextcell

  use modtrack_cell_splitting, only: splitcell

  use modarray, only: increase_array

  implicit none

  private

  ! main
  public dotracking, fillparentarr, findparents

  ! modstatistics
  public nextcell, firstcell

  contains

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

    integer(kind=2), dimension(:,:), allocatable :: cellloc

    integer :: i, j, t
    write (*,*) '.. entering tracking'

    allocate(cellloc(3,ceiling(min(0.3*(huge(1)-2),0.5*real(nx)*real(ny)*real(nt-tstart)))))
    print *, "Allocating array for storing cell locations"
print *, 'cellloc',shape(cellloc),0.3*huge(1), 0.5*real(nx)*real(ny)*real(nt-tstart)
    nullify(cell)
    do t = tstart, nt
      if(mod(t,10)==0) write (*,'(A,I10,A,I10)') "Time =", t,"  Nr of Cells", ncells
      do  j = 1, ny
        do i = 1, nx
          if (obj_mask(i,j,t)==0) then
             !write (*,*) 'Create a new cell'
            call createcell(cell)
            call newelement(i, j, t, cell, obj_mask, var_base, var_top, cellloc=cellloc)
            if (cell%nelements >= nmincells) then
              if(cell%nelements> 10000000) write (*,*) '..finalizing cell'
              call finalizecell(cell, ncells, parentarr, obj_mask, var_base, var_top, &
                                var_value, cellloc=cellloc)
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
  !! - obj_mask is set to `-1` for all datapoints for the elements in the current
  !! cell
  !!
  !! @TODO what does `ivalue` mean here?
  subroutine finalizecell(cell, ncells, parentarr, obj_mask, var_base, var_top, var_value, cellloc)
    type(celltype), pointer, intent(inout)                   :: cell
    integer, intent(out)                             :: ncells
    type(cellptr), allocatable, dimension(:,:,:), intent(inout), optional :: parentarr
    integer :: n
    integer(kind=4), dimension(:,:,:), intent(inout) :: obj_mask
    integer(kind=2), dimension(:,:,:), intent(in) :: var_base
    integer(kind=2), dimension(:,:,:), intent(in) :: var_top
    integer(kind=2), dimension(:,:,:), intent(in) :: var_value
    integer(kind=2), dimension(:,:), intent(in) :: cellloc

    if (present(parentarr)) then
      call splitcell(cell, ncells, parentarr, obj_mask, var_base, var_top, var_value, &
                     cellloc=cellloc)
    else
      allocate(cell%loc(3,cell%nelements))
      allocate(cell%value(3,cell%nelements))
      cell%loc(:,1:cell%nelements) = cellloc(:,1:cell%nelements)
      do n = 1, cell%nelements
        cell%value(ibase, n) = var_base(cellloc(1,n), cellloc(2,n), cellloc(3,n))
        cell%value(itop, n)  = var_top(cellloc(1,n), cellloc(2,n), cellloc(3,n))
        cell%value(ibase, n) = var_value(cellloc(1,n), cellloc(2,n), cellloc(3,n))
        obj_mask(cellloc(1,n), cellloc(2,n), cellloc(3,n)) = -1
      end do
    end if
    ncells = ncells + 1
    cell%id = ncells
  end subroutine finalizecell

  !> Given the current location (i,j,t) in space and time look at neighbouring
  !> points in space and time and if they satisfy the constraints from being part of
  !> the same cell:
  !!
  !! set obj_mask=-2 so that this data-point is not considered twice for multiple
  !! cells, store the position in space and time into the `cellloc` array and
  !! increment the `nelements` counter on the provided cell
  !!
  !! TODO: Looking west/east/north/south and truncating the indexing near the
  !! edge is unnecessarily costly, the current element will be checked twice
  recursive subroutine newelement(i, j, t, cell, obj_mask, var_base, var_top, cellloc)
    integer, intent(in)                                      :: i, j, t
    type(celltype),pointer, intent(inout)                    :: cell
    integer :: ii, jj, tt
    integer(kind=4), dimension(:,:,:), intent(inout) :: obj_mask
    integer(kind=2), dimension(:,:), intent(inout) :: cellloc

    integer(kind=2), dimension(:,:,:), intent(in) :: var_base
    integer(kind=2), dimension(:,:,:), intent(in) :: var_top

    cell%nelements = cell%nelements + 1
    if (mod(cell%nelements,10000000) == 0) write(*,*) 'Cell element ', cell%nelements
    cellloc(1,cell%nelements) = i
    cellloc(2,cell%nelements) = j
    cellloc(3,cell%nelements) = t

    obj_mask(i, j, t) = -2

    !Look west
    ii = i - 1
    jj = j
    tt = t
    if (ii.le.0) ii = nx
    if (obj_mask(ii,jj,tt)==0) then
      if (var_base(i,j,t) <= var_top(ii,jj,tt) .and. var_top(i,j,t) >= var_base(ii,jj,tt)) then
        call newelement(ii,jj,tt,cell, obj_mask, var_base, var_top, cellloc=cellloc)
      end if
    end if
    !Look east
    ii = i + 1
    jj = j
    tt = t
    if (ii.gt.nx) ii = 1
    if (obj_mask(ii,jj,tt)==0) then
      if (var_base(i,j,t) <= var_top(ii,jj,tt) .and. var_top(i,j,t) >= var_base(ii,jj,tt)) then
        call newelement(ii,jj,tt,cell, obj_mask, var_base, var_top, cellloc=cellloc)
      end if
    end if

    !Look north
    ii = i
    jj = j - 1
    tt = t
    if (jj.le.0) jj = ny
    if (obj_mask(ii,jj,tt)==0) then
      if (var_base(i,j,t) <= var_top(ii,jj,tt) .and. var_top(i,j,t) >= var_base(ii,jj,tt)) then
        call newelement(ii,jj,tt,cell, obj_mask, var_base, var_top, cellloc=cellloc)
      end if
    end if

    !Look south
    ii = i
    jj = j + 1
    tt = t
    if (jj.gt.ny) jj = 1
    if (obj_mask(ii,jj,tt)==0) then
      if (var_base(i,j,t) <= var_top(ii,jj,tt) .and. var_top(i,j,t) >= var_base(ii,jj,tt)) then
        call newelement(ii,jj,tt,cell, obj_mask, var_base, var_top, cellloc=cellloc)
      end if
    end if

    !Look forward
    ii = i
    jj = j
    tt =t + 1
    if (tt <= nt) then
      if (obj_mask(ii,jj,tt)==0) then
        if (var_base(i,j,t) <= var_top(ii,jj,tt) .and. var_top(i,j,t) >= var_base(ii,jj,tt)) then
          call newelement(ii,jj,tt,cell, obj_mask, var_base, var_top, cellloc=cellloc)
        end if
      end if
    end if
    !Look backward
    ii = i
    jj = j
    tt = t - 1
    if (tt >= tstart) then
      if (obj_mask(ii,jj,tt)==0) then
        if (var_base(i,j,t) <= var_top(ii,jj,tt) .and. var_top(i,j,t) >= var_base(ii,jj,tt)) then
          call newelement(ii,jj,tt,cell, obj_mask, var_base, var_top, cellloc=cellloc)
        end if
      end if
    end if
  end subroutine newelement

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
    integer(kind=2), allocatable,dimension(:), intent(in)      :: minbase
    type(cellptr), allocatable, dimension(:,:,:), intent(inout)   :: parentarr
    integer(kind=2), allocatable,dimension(:,:,:), intent(out), optional  :: base, top

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
      if (cell%nelements > n_minparentel) then
        print *, minval(cell%value(ibase,1:cell%nelements)), minbase(cell%loc(3,1))

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
    integer(kind=4), dimension(:,:,:), intent(inout) :: obj_mask

    integer :: t, count1, count2
    real    :: hlp

      !..check for empty frames
    write(0,*) 'Check frames'
    count2 = 0
    do t = 2, nt
      count1 = count(obj_mask(:,:,t)==0)
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
