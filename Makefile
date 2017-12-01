NETCDF_F_PATH=/sw/squeeze-x64/netcdf_fortran-4.2.0-static-gcc47

NF_CONFIG="$(NETCDF_F_PATH)/bin/nf-config"

FFLAGS=$(shell $(NF_CONFIG) --fflags --flibs --cflags) 
FC=$(shell $(NF_CONFIG) --fc)
NETCDF_INCDIR=$(shell $(NF_CONFIG) --prefix)/include


all: cloud_tracking

modnetcdf: modnetcdf.f90
	$(FC) $(FFLAGS) -c $(addsuffix .f90,$@)

cloud_tracking: cloud_tracking.f90 modnetcdf
	$(FC) $(addsuffix .f90,$@) $(FFLAGS) modnetcdf.o -o bin/$@
