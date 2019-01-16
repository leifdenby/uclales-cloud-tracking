module modstatistics
 use tracking_common, only: celltype

 use tracking_common, only: dt, dx, dy
 use tracking_common, only: nx, ny
 use tracking_common, only: ibase, itop, ivalue
 use tracking_common, only: tstart, tend

 use modtrack, only: nextcell, firstcell

 implicit none

 contains
  subroutine find_max_number_of_relatives(cell, n)
    use modtrack, only: nextcell, firstcell

    type(celltype), pointer, intent(inout) :: cell
    integer :: n
    integer :: iret

    iret = firstcell(cell)
    do
      if (iret == -1) exit
      n = max(n, max(cell%nparents, cell%nchildren))
      iret = nextcell(cell)
    end do
  end subroutine

  subroutine dostatistics(cell, ncells, icell, fid, ivar, heightzero, heightrange, distmax, valzero, valrange, ovarstem)
    use modnetcdf, only: netcdfvar, fillvalue_r, fillvalue_i
    use modnetcdf, only: xnc, ync, tnc
    use modnetcdf, only: define_ncvar, write_ncvar, define_ncdim
    use netcdf, only: nf90_float, nf90_int
    use constants, only: real_maxval

    type(celltype), pointer, intent(inout)    :: cell
    integer, intent(in)               :: ncells
    integer, intent(in)                       :: fid, icell
    real, intent(in)                          :: heightzero, heightrange, distmax,valzero, valrange
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
    integer, dimension(nx, ny, tstart:tend) :: slab
    real, dimension(:,:), allocatable :: base, top, area, vol, val, xcenter, ycenter, maxarea, maxarealoc
    !real, dimension(:,:), allocatable :: recon
    real, dimension(:), allocatable   :: duration, mintime, maxtime
    real :: dz = 25, rbase, rtop
    integer, dimension(:,:), allocatable :: relatives
    integer, dimension(:), allocatable   :: id, nrelatives
    real :: time

    integer :: max_num_relatives = 0

    write (*,*) '.. entering statistics for ', ncells, " ", trim(ovarstem%name), 's'

    !Loop over the cells  - fill the cell-length distribution and the xyt slab
    slab(1:nx, 1:ny, tstart:tend) = 0
    slab(1:nx, 1:ny, tstart:tend) = -1000000000
    slab(1:nx, 1:ny, tstart:tend) = fillvalue_i
    bucketsize = 0
    bucket_max = 0
    bucket_min = 0
    bucketname = (/'sm','la','hu'/)
    bucketlongname = (/'Small','Large','Huge '/)

    call find_max_number_of_relatives(cell, max_num_relatives)

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
    relnc%dim(1)   = max_num_relatives
    call define_ncdim(fid, relnc, nf90_int)
    call write_ncvar(fid, relnc, (/(i, i=1,relnc%dim(1))/))

    n = 0
    allocate(tlength(ncells))
    tlength(:) = 0
    iret = firstcell(cell)
    do
      if (iret == -1) exit
      n = n + 1
      tmin = minval(cell%loc(3,1:cell%n_points))
      tmax = maxval(cell%loc(3,1:cell%n_points))
      tlength(n) = tmax - tmin + 1
      iret = nextcell(cell)
    end do
    if (n < ncells) then
       print *, "Error: there appear to be fewer cells than `ncells` suggests"
       print *, n, ncells
       call exit(1)
    else if (n > ncells) then
       print *, "Error: there appear to be more cells than `ncells` suggests"
       print *, n, ncells
       call exit(1)
    endif

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
        !allocate(recon(ceiling(distmax/dz), bucket_max(n)))
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
          tmin = minval(cell%loc(3,1:cell%n_points))
          tmax = maxval(cell%loc(3,1:cell%n_points))
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
            !recon = 0.
!             icenter = 0
!             jcenter = 0
!             npts  = 0

            do nel = 1, cell%n_points
              i = cell%loc(1,nel)
              j = cell%loc(2,nel)
              t = cell%loc(3,nel)
              tt = t - tmin + 1
              if (cell%id<0) then
                write (*,*) 'DANGER: Cell id < 0'
              else
                slab(i,j,t) = cell%id
              end if
              rbase = real(cell%value(ibase,nel)) * heightrange/real_maxval + heightzero
              rtop  = real(cell%value(itop,nel)) * heightrange/real_maxval + heightzero
              base(tt, nn) = min(base(tt, nn), rbase)
              top (tt, nn) = max(top (tt, nn), rtop)
              vol (tt, nn) = vol(tt, nn) + rtop-rbase
              area(tt, nn) = area(tt, nn) + 1.
              !recon(floor(rbase/dz)+1:floor(rtop/dz)+1,tt) = recon(floor(rbase/dz)+1:floor(rtop/dz)+1,tt) + 1.
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
              !maxarea(tt,nn) = maxval(recon(:,tt)) *dx*dy
              !maxarealoc(tt,nn) = maxloc(recon(:,tt),1)*dz
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

        deallocate(base, top, area, maxarea, maxarealoc, vol, val)
        !deallocate(recon)
        deallocate(icenter,jcenter, xcenter, ycenter,ianchor, janchor, npts)
        deallocate(duration, mintime, maxtime, id)



        write (*,*) '..Relationships'
        allocate(nrelatives(bucketsize(n)))
        allocate(relatives(max_num_relatives,bucketsize(n)))
        relatives = fillvalue_i

        print *, "nrel", shape(relatives)


        nn   = 0
        iret = firstcell(cell)
        do
          if (iret == -1) exit
          tmin = minval(cell%loc(3,1:cell%n_points))
          tmax = maxval(cell%loc(3,1:cell%n_points))
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
          tmin = minval(cell%loc(3,1:cell%n_points))
          tmax = maxval(cell%loc(3,1:cell%n_points))
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
            tmin = minval(cell%loc(3,1:cell%n_points))
            tmax = maxval(cell%loc(3,1:cell%n_points))
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
            tmin = minval(cell%loc(3,1:cell%n_points))
            tmax = maxval(cell%loc(3,1:cell%n_points))
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
          tmin = minval(cell%loc(3,1:cell%n_points))
          tmax = maxval(cell%loc(3,1:cell%n_points))
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
end module modstatistics
