#! /bin/sh -v

# Define the fortran compiler and options

#include ./Makefile.conf

FC=ifort
FLAGS=-convert big_endian

# Set these to match your system's installation paths
NCDIAG_DIR =
NETCDF_DIR  =

INCLD = -I$(NCDIAG_DIR)/include -I$(NETCDF_DIR)/include
LIBS  = -L$(NCDIAG_DIR)/lib -lncdiag \
        -L$(NETCDF_DIR)/lib -lnetcdff -lnetcdf

all: kinds read_diag_conv_rtma

read_diag_conv_rtma: read_diag_conv_rtma.o kinds.o
	${FC} -o read_diag_conv_rtma.exe ${FLAGS} read_diag_conv_rtma.o kinds.o ${LIBS}
	rm read_diag_conv_rtma.o

kinds: kinds.F90
	$(FC) $(FLAGS) -c kinds.F90

.SUFFIXES : .f90  .o

.f90.o :
	${FC} ${FLAGS} ${INCLD} -c $<

clean:
	/bin/rm -f *.o *.exe *.mod
