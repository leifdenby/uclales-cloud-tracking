1. What input does the cloud tracking code expect?

    ./tracking {var_name} {t_start} {nt} {criterion}

    var_name: variable name
    t_start: starting timestep
    nt: number of timesteps (I think)
    criterion: what to search for
        - 'core': need `lwp` and `core`
        - 'cloud': need `lwp`
        - 'liquid': need `lwp` and `rain`
        - `thermal`: need `trcpath`
        - `rain`: need `rain`
        - `all`: need all above variables

    Expects data to be in files with all horizontally projected data for
    multiple timesteps for single variable `var_name`:

    `{experiment_name}.out.xy.{var_name}.nc`

    - How do I create this file?

2. How does it load the input and into what?

- All data is loaded into `var(nx, ny, tstart:nt, nvar)` where `nvar` is the
index that has been assigned to the currently loaded variable, the value of
`nvar` for each loaded variable is defined in the argument parsing routine.
`tstart:nt` indicates that the indexing in the allocated array runs from
`tstart` to `nt`

- All data is offset and rescaled using an offset and range value, given as
  `{var_name}zero` and `{var_name}range` for variable `var_name`.
  `{var_name}range` is given by partitioning the maximum int value through the
  absolute range of values specified for `{var_name}` with `{var_name}max` and
  `{var_name}min`, and `{var_name}zero` is set as the halfway-point between
  `{var_name}min` and `{var_name}max`. *NB*: this rescaling means that the
  effective resolution of the output is dependent on choosing `{var_name}min`
  and `{var_name}max` carefully.

TODO: `{var_name}min` and `{var_name}max` could be determined from the data I think.

- data is loaded for the chosen timesteps range using the `nf90_get_var`
  routine of the netcdf library

3. Where does the analysis start?

## analysis of `core`s

    - starts off as `-1`

- during analysis it appears that "cells" are stored in a linked list, starting
  with `thermal` which is a pointer to a object of type `celltype`


## `clouds` analysis

`if (lcloud) then`, l. 2001 in `tracking_time.f90`

## global variables

### `cellloc`

- in `newelement` (i,j,t) (indexes in space and time) and assigned to
  `newelement(1,...)`, `newelement(2,...)` and `newelement(3,...)` respectively

### `bool`

- allocated to have shape `bool(nx, ny, tstart:nt)`
- initiated to -1
- in "core" analysis:
    - cells where `lwp` > `i_lwpthres` => initiated to `bool` = 0
- in `newelement` subroutine:
    - set to `-2` for `(i,j,t)`-indexes passed in
    - required to be == 0 for a cell to be consider for a "new element"
- appears to indicate type of a particular cell

## functions and subroutines

- firstcell ({->celltype})
- nextcell ({->celltype})
- check
- checkframes
    - Appears to check for non-consistent changes in `bool`s values in
      successive timesteps. `bool==0` is appearently expected not to drop to 0
      if it was > 100, and should also not change by more than 75%
- createcell ({->celltype})
    - if passed an associated `celltype` pointer a "neighbouring" (assigned to
      `next`) `celltype` instance is created, with same `head` and `next` as
      the provided `celltype` instance
    - if not associated and new `celltype` instance (with head pointing to
      itself) is created and the provided pointer is associated to it
- delete_all ({->celltype})
    - use `deletecell` to remove all elements in `next`/`previous` linked-list
      that the provided `celltype` instance exist in
- deletecell ({->celltype})
    - delete single element in `celltype` `next`/`previous` linked-list,
      deleting also all allocated "parents", "children", "splitters" celltypes
- dostatistics
- dotracking
- fillparentarr
  !> Create a mapping from all datapoints in space and time to the cell
  !  associated with that datapoint for all cells where the minimum "cloud-base"
  !  value is less than minimum base height provided (this is actually and array
  !  that spans all timesteps, but only the first time-step index of first
  !  element in the cell is used).
- finalizecell
- findneighbour
- findnrsplitters
- findparents
- newelement
  !> Given the current location (i,j,t) in space and time look at neighbouring
  ! points in space and time and if they satisfy the constraints from being part of
  ! the same cell:
  ! - 
  ! set bool=-2 so that this data-point is not considered twice for multiple
  ! cells, store the position in space and time into the `cellloc` array and
  ! increment the `nelements` counter on the provided cell
- newpassive

- splitcell

*dynamic array allocation*

- increase_array_i
- increase_array_p
- increase_array_r

*IO*

- define_ncdim
- define_ncvar
- inquire_ncvar
- read_ncvar_1D
- read_ncvar_2D
- read_ncvar_3D
- write_ncvar_1D_i
- write_ncvar_1D_i64
- write_ncvar_1D_r
- write_ncvar_2D_i
- write_ncvar_2D_i64
- write_ncvar_2D_r
- write_ncvar_3D_i
- write_ncvar_3D_i64
- write_ncvar_3D_r

## Datatypes

`celltype`:
- id
- nelements
- loc
- value
    - allocated and set in `finalizecell`, `value(3, cell%nelements)`
        - value(1, element_index) = variable value for cloud-base at element_index (i,j,k)
        - value(2, element_index) = variable value for cloud-top at element_index (i,j,k)
        - value(3, element_index) = variable value for "???" at element_index (i,j,k)

- nsplitters (int)
- splitters ([celltype])

- nparents (int)
- parents ([celltype])

- nparents (int)
- parents ([celltype])

- nchildren (int)
- children ([celltype])

- cloudsystemnr
- cloudtype
    - starts of as `-1` in `createcell` subroutine

- next (-> celltype)
- previous (-> celltype)
- head (-> celltype)


# Questions

- what is the `nchunck` during loading?
- what is `minparental`?
- `nmincells` = 4, argument to `dotracking`?
- whats with the 10000000 the crops up everywhere?
  - in `dotracking` `finalizecell` is call if `cell%nelements > 1000000`

# Ideas

- move:
    - the creation of a mapping from data-points in space and time
    - from: a subroutine called within main
    - to: a subroutine contained somewhere else and called from the
      `dotracking` routine if needed (i.e. if a "parent" set of cells is passed
      in)
