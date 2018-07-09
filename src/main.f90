!> Main entrypoint of cloud tracking code
!>
program tracking
  use modnetcdf, only: netcdfvar, fillvalue_i16, check
  use modnetcdf, only: read_ncvar, inquire_ncvar, define_ncdim, write_ncvar
  use modnetcdf, only: tnc, xnc, ync

  use netcdf, only: nf90_close, nf90_global, nf90_nowrite
  use netcdf, only: nf90_hdf5
  use netcdf, only: nf90_put_att, nf90_create, nf90_enddef, nf90_open

  use tracking_data, only: tracked_cores, tracked_clouds, tracked_rain, tracked_thermals
  use tracking_data, only: icore, irain, ilwp, ithermal
  use tracking_data, only: ncores, nrains, nthermals, nclouds
  use tracking_data, only: nvar

  use tracking_common, only: cellptr
  use tracking_common, only: dt, dx, dy
  use tracking_common, only: nt, nx, ny
  use tracking_common, only: tstart, simulation_id
  use tracking_common, only: INSIDE_OBJECTS, OUTSIDE_OBJECTS

  use modtrack, only: dotracking, findparents
  use modtrack, only: fillparentarr

  use field_loader, only: read_named_input
  use field_loader, only: var
  use field_loader, only: var_ibase, var_itop

  use constants, only: nchunk
  use constants, only: nmincells, nmincells_cloud
  use constants, only: i_corethres, i_lwpthres, i_thermthres, i_rwpthres
  use constants, only: distrange, distzero
  use constants, only: lwpzero, rwpzero, corezero, thermzero
  use constants, only: lwpthres, rwpthres, corethres, thermthres
  use constants, only: lwprange, rwprange, corerange, thermrange
  use constants, only: maxheight

  use modstatistics, only: dostatistics

  implicit none

  ! runtime variables
  type(cellptr), dimension(:,:,:), allocatable :: parentarr

  logical :: lcore = .false., lcloud = .false., lthermal = .false., lrain = .false., lrwp = .false.
  integer :: n
  integer :: fid
  type(netcdfvar), dimension(10) :: ivar
  type(netcdfvar) :: ovar
  integer(kind=2), dimension(:,:,:), allocatable :: var_base, var_top
  integer(kind=2), dimension(:), allocatable :: minbasecloud, minbasetherm

  !> stores object mask and state of processing for a given point, as well as later the unique index of a cell's index that splits
  !> the cell a given point belongs to
  integer(kind=4), dimension(:,:,:), allocatable :: obj_mask

  integer :: var_ivalue  ! indicates which field of `var` we are currently using

  if (command_argument_count() == 0) then
     print *, "./tracking [filename-base] [starting timestep] [final timestep] [analysis variables]"
     print *, ""
     print *, "analysis variables: `core`, `cloud`, `liquid`??, `thermal`, `rain` and `all`"
     print *, "analysis variables defines which fields will be analysed"
     call exit(-1)
  endif

  call parse_command_arguments()
  call summarise_active_fields()
  call check_for_required_files(simulation_id)
  call setup_output(simulation_id)
  call read_input(simulation_id)

  print *, "Allocating 3D arrays for storage"

  print *, "var =>", nx, ny, tstart, nt, nvar
  allocate(var(nx, ny, tstart:nt, nvar))
  allocate(var_base(nx, ny, tstart:nt))
  allocate(var_top(nx, ny, tstart:nt))

  print *, "parentarr =>", nx, ny, tstart, nt
  allocate(parentarr(nx, ny, tstart:nt))
  allocate(obj_mask(nx, ny, tstart:nt))

  ! Identify and track thermals using Couvreux "radioactive scalar", from which the min-height,
  ! max-height and the column-integrated value of this scalar
  if (lthermal) then
    write (*,*) 'Thermals....'
    call read_named_input(var(:,:,:,ithermal), ivar(ithermal))
    var_ivalue = ithermal
    print *, ">>", nvar, ithermal

    ivar(var_ibase)%name  = 'trcbase'
    ivar(var_itop)%name  = 'trctop'
    call read_named_input(var_base, ivar(var_ibase))
    call read_named_input(var_top, ivar(var_itop))

    allocate(minbasetherm(tstart:nt))
    minbasetherm = (100.-distzero)/distrange

    obj_mask = OUTSIDE_OBJECTS
    where (var(:,:,:,var_ivalue) > i_thermthres)
      obj_mask = INSIDE_OBJECTS
    end where
    call dotracking(tracked_thermals, nthermals,nmincells, obj_mask, var_base, var_top, var(:,:,:,var_ivalue))
  end if


  ! track cloud-cores and clouds
  if (lcore .or. lcloud) then
    ! load lwp + rwp => var(:,:,:,var_ivalue), rescaling by lwpzero and lwprange
    if (.not. lrwp) then
      call read_named_input(var(:,:,:,ilwp), ivar(ilwp))
    else
      call read_named_input(var(:,:,:,ilwp), ivar(ilwp), ivar(irain))
    endif

    ! load cldbase => var(:,:,:,var_ibase)
    ! load cldtop  => var(:,:,:,var_itop)
    ivar(var_ibase)%name  = 'cldbase'
    ivar(var_itop)%name  = 'cldtop'
    call read_named_input(var_base, ivar(var_ibase))
    call read_named_input(var_top, ivar(var_itop))

    ! track the cloud cores (regions of buoyancy) using a critical value for buoyancy as the mask
    ! and the cloud-top and cloud-base height to decide connectivity
    if (lcore) then
      write (*,*) 'Cores....'
      var_ivalue = 4

      call read_named_input(var(:,:,:,var_ivalue), ivar(icore))

      obj_mask = OUTSIDE_OBJECTS
      allocate(minbasecloud(tstart:nt))
      minbasecloud = 0

      do n=tstart, nt
        ! define minimum cloud-base height as
        !   minbasecloud = z_top_max + z_base_min_half
        ! where
        !   z_top_max : rescaled domain-wide maximum cloud-top height
        !   z_base_min_half : half of rescaled domain-wide minimum cloud-base height (where the cloud-base heigth is defined)
        ! (20/06/2018, Leif: it's unclear to me why using z_top_max isn't enough)
        minbasecloud(n) = (maxval(var_top(:,:,n)) + minval(var_base(:,:,n),var_base(:,:,n)>fillvalue_i16) )/2

        if (any(var(:,:,n,3) > i_lwpthres)) then
          ! where
          !       LWP is above threshold
          !   and buoyancy is above threshold
          !   and cloud-base height is below cloud-top height
          where (var(:,:,n,3) > i_lwpthres &
              .and. var(:,:,n,var_ivalue) > i_corethres &
              .and. var_base(:,:,n) < minbasecloud(n))
            obj_mask(:,:,n) = INSIDE_OBJECTS
          end where
        end if
      end do

      !       call checkframes
      call dotracking(tracked_cores, ncores,nmincells, obj_mask, var_base, var_top, var(:,:,:,var_ivalue))
    end if
    ! lwp field is the 3 entry in the data array
    var_ivalue = 3

    ! Track cloud using a critical value of column-integrated liquid as the mask and cloud-top and
    ! cloud-base height as the mask. If cloud-cores have been tracked too the clouds are split by
    ! their cores
    if (lcloud) then
      write (*,*) 'Clouds....'
      obj_mask = OUTSIDE_OBJECTS
      where (var(:,:,:,ilwp) > i_lwpthres)
        obj_mask = INSIDE_OBJECTS
      end where

      if (lcore) then
        call fillparentarr(tracked_cores, minbasecloud, parentarr)
        call dotracking(tracked_clouds, nclouds,nmincells_cloud, obj_mask, var_base, var_top, var(:,:,:,var_ivalue), parentarr)

        !call fillparentarr(tracked_clouds, minbasecloud, parentarr)
        !call findparents(tracked_cores, parentarr, var_base, var_top)
      else
        call dotracking(tracked_clouds, nclouds,nmincells_cloud, obj_mask, var_base, var_top, var(:,:,:,var_ivalue))
      end if
      if (lthermal) then
        ! NB: `var_base` and `var_top` are overwritten here, we reuse them to save memory
        call fillparentarr(tracked_thermals, minbasetherm, parentarr, var_base, var_top)
        call findparents(tracked_clouds, parentarr, var_base, var_top)
        deallocate(var_base,var_top)
      end if
    end if
  end if

  if (lrain) then
    write (*,*) 'Rain....'
    call read_named_input(var(:,:,:,irain), ivar(irain))

    ivar(var_ibase)%name  = 'rwpbase'
    ivar(var_itop)%name  = 'rwptop'
    call read_named_input(var_base, ivar(var_ibase))
    call read_named_input(var_top, ivar(var_itop))

    !Loop over rain patches
    obj_mask = OUTSIDE_OBJECTS
    where (var(:,:,:,var_ivalue) > i_rwpthres)
      obj_mask = INSIDE_OBJECTS
    end where

    call dotracking(tracked_rain, nrains,nmincells, obj_mask, var_base, var_top, var(:,:,:,var_ivalue))
    if (lcloud) then
      if (.not. allocated(minbasecloud)) allocate(minbasecloud(tstart:nt))
      minbasecloud = huge(1_2)
      ! NB: `var_base` and `var_top` are overwritten here, we reuse them to save memory
      call fillparentarr(tracked_clouds, minbasecloud, parentarr, var_base, var_top)
      call findparents(tracked_rain, parentarr, var_base, var_top)
    end if
  end if

  if (lthermal) then
    ovar%name     = 'thrm'
    ovar%longname = 'Thermal'
    call dostatistics(tracked_thermals, nthermals, ithermal, fid, ivar, distzero, distrange, maxheight,thermzero, thermrange, ovar)
!     call delete_all(tracked_thermals)
  end if

  if (lcore) then
    ovar%name     = 'core'
    ovar%longname = 'Core'
    call dostatistics(tracked_cores, ncores, icore, fid, ivar, distzero, distrange, maxheight,corezero, corerange, ovar)
!     call delete_all(core)
  end if

  if (lcloud) then
    ovar%name     = 'cloud'
    ovar%longname = 'Cloud'
    call dostatistics(tracked_clouds, nclouds, ilwp, fid, ivar, distzero, distrange, maxheight,lwpzero, lwprange, ovar)
!     call delete_all(tracked_clouds)
  end if

  if (lrain) then
    ovar%name     = 'rain'
    ovar%longname = 'Rain patch'
    call dostatistics(tracked_rain, nrains, irain, fid, ivar, distzero, distrange, maxheight, rwpzero, rwprange, ovar)
!     call delete_all(tracked_rain)
  end if
  call check ( nf90_close(fid))
  write (*,*) '..Done'

  contains
     subroutine check_file_exists(filename)
        character(len=100), intent(in) :: filename
        logical :: file_exists
        inquire(file=trim(filename), exist=file_exists)

        if (.not. file_exists) then
           print *, trim(filename)//" is missing"
           call exit(1)
        endif
     end subroutine check_file_exists

     subroutine check_for_required_files(fsimulation_id)
        character(len=100), intent(in) :: fsimulation_id
        character(100) :: fname

        if (icore /= 0) then
          fname = trim(fsimulation_id)//'core.nc'
          call check_file_exists(fname)
        endif

        if (ilwp /= 0) then
          fname = trim(fsimulation_id)//'lwp.nc'
          call check_file_exists(fname)
        endif
     end subroutine check_for_required_files

     subroutine setup_output(simulation_id)
        character(len=100), intent(in) :: simulation_id

        call check ( nf90_create(path=trim(simulation_id)//'track.nc', cmode=NF90_HDF5, ncid=fid))
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

     subroutine read_input(simulation_id)
        character(len=100), intent(in) :: simulation_id
        character(100) :: filename
        real, allocatable, dimension(:) :: x, y, t
        integer :: finput

        filename = trim(simulation_id)//trim(ivar(nvar)%name)//'.nc'
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
        print *, nt, "timesteps in file"
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
        call check(nf90_close(finput))
        write (*,*) nx, ny, nt, tstart,nvar
     end subroutine read_input

     ! Parse command line arguments
     ! `./tracking_time [filename-base] [starting timestep] [final timestep] [analysis variables]
     !
     ! analysis variables: `core`, `cloud`, `liquid`??, `thermal`, `rain` and `all`
     ! analysis variables defines which fields will be analysed
     subroutine parse_command_arguments()
        character(100) :: criterion, ctmp
        integer :: ub, lb, vlength

        call get_command_argument(1,simulation_id)
        simulation_id = trim(simulation_id)//'.out.xy.'
        call get_command_argument(2,ctmp)
        read(ctmp,'(i4)') tstart
        call get_command_argument(3,ctmp)
        read(ctmp,'(i4)') nt
        call get_command_argument(4,criterion)

        nvar = 2
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
        case default
           print *, trim(adjustl(criterion(lb:ub)))//" isn't a valid option"
           call exit(2)
        end select
        end do
        if (ilwp     > 0) ivar(ilwp)%name     = "lwp"
        if (icore    > 0) ivar(icore)%name    = 'core'
        if (ithermal > 0) ivar(ithermal)%name = 'trcpath'
        if (irain    > 0) ivar(irain)%name    = 'rwp'
        ! End parse command line arguments

        print *, "going to use", nvar, "variables"
     end subroutine parse_command_arguments

     subroutine summarise_active_fields
        integer :: n
        character(len=10) :: n_str
        print *, "Active fields:"

        do n=1, nvar
          write (n_str, "(I1)") n
          print *, "("//trim(n_str)//") => "//trim(ivar(n)%name)
        enddo

     end subroutine summarise_active_fields

end program tracking
