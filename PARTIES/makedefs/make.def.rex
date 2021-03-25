# -*- mode: makefile -*-
# Machine: rex

# Compiler executable
compiler_cmd := mpicc
# Compiler flags
#   cflags - flags to give the compiler
#   optmz  - optimization (default)
#   debug  - debugging
#   libs   - libraries needed to be linked to executable
cflags      +=
optmz_flags := -O3
debug_flags := -g
libs        += -lm -ldl -lfftw3 -lhdf5

# Location of packages
fftw_dir := /opt/fftw-3.3.6-pl1
hdf5_dir := /opt/hdf5-1.8.18
mpi_dir  := #/opt/openmpi-2.0.2

# FFTW
ifneq ($(strip $(fftw_dir)),) # If not empty
  lflags   += -L$(fftw_dir)/lib
  libs     += -Wl,-rpath=$(fftw_dir)/lib
  includes += -I$(fftw_dir)/include
endif

# HDF5
ifneq ($(strip $(hdf5_dir)),) # If not empty
  lflags   += -L$(hdf5_dir)/lib
  libs     += -Wl,-rpath=$(hdf5_dir)/lib
  includes += -I$(hdf5_dir)/include
endif

# MPI
ifneq ($(strip $(mpi_dir)),) # If not empty
  lflags   += -L$(mpi_dir)/lib/openmpi
  libs     += -Wl,-rpath=$(mpi_dir)/lib/openmpi
  includes += -I$(mpi_dir)/include
  compiler := $(mpi_dir)/bin/$(compiler_cmd)
else
  compiler := $(compiler_cmd)
endif
