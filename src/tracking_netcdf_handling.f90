module modnetcdf
  use netcdf
  implicit none
  real(kind=4),    parameter :: fillvalue_r = -1e9
  integer, parameter :: fillvalue_i = -1000000000
  integer(kind=8), parameter :: fillvalue_i64 = -1000000000
  integer(kind=2), parameter :: fillvalue_i16 = -huge(1_2)
  integer :: deflate_level = 1
  type netcdfvar
     integer       :: id
     integer, dimension(:), allocatable :: dim, dimids
     character(80) :: name,longname,units
  end type netcdfvar

  type(netcdfvar) :: xnc, ync, tnc

  interface write_ncvar
    module procedure write_ncvar_1D_r
    module procedure write_ncvar_2D_r
    module procedure write_ncvar_3D_r
    module procedure write_ncvar_1D_i
    module procedure write_ncvar_2D_i
    module procedure write_ncvar_3D_i
    module procedure write_ncvar_1D_i64
    module procedure write_ncvar_2D_i64
    module procedure write_ncvar_3D_i64
  end interface write_ncvar

  interface read_ncvar
    module procedure read_ncvar_1D
    module procedure read_ncvar_2D
    module procedure read_ncvar_3D
  end interface read_ncvar
contains

  subroutine inquire_ncvar(fid, ncvar)
    type(netcdfvar) :: ncvar
    integer :: fid, n, ndims , iret
    call check ( nf90_inq_varid(fid,trim(ncvar%name),ncvar%id) )
    call check ( nf90_inquire_variable(fid,ncvar%id, ndims=ndims  ))
    if (allocated(ncvar%dim   )) deallocate(ncvar%dim)
    if (allocated(ncvar%dimids)) deallocate(ncvar%dimids)
    allocate(ncvar%dim(ndims))
    allocate(ncvar%dimids(ndims))
    call check ( nf90_inquire_variable(fid,ncvar%id, dimids = ncvar%dimids))
    do n = 1, ndims
      call check (nf90_inquire_dimension(fid, ncvar%dimids(n), len = ncvar%dim(n)))
    end do
    iret = nf90_get_att(fid, ncvar%id, 'longname',ncvar%longname)
    if (iret /= nf90_noerr) then
      ncvar%longname = ''
    end if
    iret = nf90_get_att(fid, ncvar%id, 'units',ncvar%units)
    if (iret /= nf90_noerr) then
      ncvar%units = ''
    end if

  end subroutine inquire_ncvar

  subroutine read_ncvar_1D(fid,ncvar,array,start,count)
    type(netcdfvar) :: ncvar
    real, dimension(:) :: array
    integer,dimension(1), optional, intent(in) :: start,count
    integer         :: fid

    !..open file and read data
    call check ( nf90_inq_varid(fid,trim(ncvar%name),ncvar%id) )
    if (present(start)) then
      call check ( nf90_get_var(fid,ncvar%id,array,start = start, count = count) )
    else
      call check ( nf90_get_var(fid,ncvar%id,array) )
    end if

  end subroutine read_ncvar_1D

  subroutine read_ncvar_2D(fid,ncvar,array,start,count)
    type(netcdfvar) :: ncvar
    real, dimension(:,:) :: array
    integer,dimension(2), optional, intent(in) :: start,count
    integer         :: fid

    !..open file and read data
    call check ( nf90_inq_varid(fid,trim(ncvar%name),ncvar%id) )
    if (present(start)) then
      call check ( nf90_get_var(fid,ncvar%id,array,start = start, count = count) )
    else
      call check ( nf90_get_var(fid,ncvar%id,array) )
    end if

  end subroutine read_ncvar_2D

  subroutine read_ncvar_3D(fid,ncvar,array,start,count)
    type(netcdfvar) :: ncvar
    real, dimension(:,:,:) :: array
    integer,dimension(3), optional, intent(in) :: start,count
    integer         :: fid

    !..open file and read data
    call check ( nf90_inq_varid(fid,trim(ncvar%name),ncvar%id) )
    if (present(start)) then
      call check ( nf90_get_var(fid,ncvar%id,array,start = start, count = count) )
    else
      call check ( nf90_get_var(fid,ncvar%id,array) )
    end if
  end subroutine read_ncvar_3D

  subroutine write_ncvar_1D_r(fid,ncvar,array)
    type(netcdfvar) :: ncvar
    real, dimension(:) :: array
    integer         :: fid

    !..open file and read data
    call check ( nf90_inq_varid(fid,trim(ncvar%name),ncvar%id) )
    call check ( nf90_put_var(fid,ncvar%id,array) )
  end subroutine write_ncvar_1D_r

  subroutine write_ncvar_2D_r(fid,ncvar,array)
    type(netcdfvar) :: ncvar
    real, dimension(:,:) :: array
    integer         :: fid

    !..open file and read data
    call check ( nf90_inq_varid(fid,trim(ncvar%name),ncvar%id) )
    call check ( nf90_put_var(fid,ncvar%id,real(array)) )
  end subroutine write_ncvar_2D_r

  subroutine write_ncvar_3D_r(fid,ncvar,array)
    type(netcdfvar) :: ncvar
    real, dimension(:,:,:) :: array
    integer         :: fid

    !..open file and read data
    call check ( nf90_inq_varid(fid,trim(ncvar%name),ncvar%id) )
    call check ( nf90_put_var(fid,ncvar%id,array) )
  end subroutine write_ncvar_3D_r

  subroutine write_ncvar_1D_i(fid,ncvar,array)
    type(netcdfvar) :: ncvar
    integer, dimension(:) :: array
    integer         :: fid

    !..open file and read data
    call check ( nf90_inq_varid(fid,trim(ncvar%name),ncvar%id) )
    call check ( nf90_put_var(fid,ncvar%id,array) )
  end subroutine write_ncvar_1D_i

  subroutine write_ncvar_2D_i(fid,ncvar,array)
    type(netcdfvar) :: ncvar
    integer, dimension(:,:) :: array
    integer         :: fid

    !..open file and read data
    call check ( nf90_inq_varid(fid,trim(ncvar%name),ncvar%id) )
    call check ( nf90_put_var(fid,ncvar%id,array) )
  end subroutine write_ncvar_2D_i

  subroutine write_ncvar_3D_i(fid,ncvar,array)
    type(netcdfvar) :: ncvar
    integer, dimension(:,:,:) :: array
    integer         :: fid

    !..open file and read data
    call check ( nf90_inq_varid(fid,trim(ncvar%name),ncvar%id) )
    call check ( nf90_put_var(fid,ncvar%id,array) )
  end subroutine write_ncvar_3D_i

  subroutine write_ncvar_1D_i64(fid,ncvar,array)
    type(netcdfvar) :: ncvar
    integer(kind=8), dimension(:) :: array
    integer         :: fid

    !..open file and read data
    call check ( nf90_inq_varid(fid,trim(ncvar%name),ncvar%id) )
    call check ( nf90_put_var(fid,ncvar%id,array) )
  end subroutine write_ncvar_1D_i64

  subroutine write_ncvar_2D_i64(fid,ncvar,array)
    type(netcdfvar) :: ncvar
    integer(kind=8), dimension(:,:) :: array
    integer         :: fid

    !..open file and read data
    call check ( nf90_inq_varid(fid,trim(ncvar%name),ncvar%id) )
    call check ( nf90_put_var(fid,ncvar%id,array) )
  end subroutine write_ncvar_2D_i64

  subroutine write_ncvar_3D_i64(fid,ncvar,array)
    type(netcdfvar) :: ncvar
    integer(kind=8), dimension(:,:,:) :: array
    integer         :: fid

    !..open file and read data
    call check ( nf90_inq_varid(fid,trim(ncvar%name),ncvar%id) )
    call check ( nf90_put_var(fid,ncvar%id,array) )
  end subroutine write_ncvar_3D_i64

  subroutine define_ncdim(fid,ncvar, dimtype)
    type(netcdfvar)   :: ncvar
    integer           :: fid, iret, xtype
    integer, optional :: dimtype

    if (present(dimtype)) then
      xtype = dimtype
    else
      xtype = nf90_float
    end if
    !..open file and read data
    call check ( nf90_redef(fid))
    iret = nf90_def_dim(fid, trim(ncvar%name), ncvar%dim(1), ncvar%dimids(1))
    if (iret == nf90_enameinuse) then
      call check ( nf90_inq_varid(fid,trim(ncvar%name),ncvar%dimids(1)) )
    else
      call check(iret)
    end if
    call check ( nf90_enddef(fid))
    call define_ncvar(fid, ncvar, xtype)

  end subroutine define_ncdim

  subroutine define_ncvar(fid,ncvar, xtype)
         INCLUDE 'netcdf.inc'
    type(netcdfvar) :: ncvar
    integer         :: fid, iret, xtype
    !..open file and read data
    call check ( nf90_redef(fid))
    iret =  nf90_def_var(fid, name = trim(ncvar%name),xtype=xtype, dimids = ncvar%dimids, &
                         varid = ncvar%id, deflate_level = deflate_level)
    if (iret == nf90_enameinuse) then
      call check ( nf90_inq_varid(fid,trim(ncvar%name),ncvar%id) )
    else
      call check(iret)
      call check ( nf90_put_att(fid, ncvar%id, 'longname', ncvar%longname))
      call check ( nf90_put_att(fid, ncvar%id, 'units', ncvar%units))
      select case(xtype)
      case(nf90_int)
        call check ( nf90_put_att(fid, ncvar%id, '_FillValue', fillvalue_i))
      case(nf90_int64)
        call check (nf_put_att_int(fid, ncvar%id, '_FillValue', nf90_int64, 1, (/fillvalue_i/)))
      case(nf90_float)
        call check ( nf90_put_att(fid, ncvar%id, '_FillValue', fillvalue_r))
      end select
    end if
    call check ( nf90_enddef(fid))

  end subroutine define_ncvar

  subroutine check(status, allowed)
    integer, intent (in) :: status
    integer, dimension(:), intent (in), optional :: allowed
    if(status /= nf90_noerr) then
      write(0,*)  status, trim(nf90_strerror(status))
      if (present(allowed)) then
        if (any(allowed == status)) return
      end if
      stop 2
    end if
  end subroutine check

end module modnetcdf
