
# Load system information
ifdef MY_MACHINE
include make.def.$(MY_MACHINE)
else
include make.def.generic
endif

# Variables from system information to be used by other makefiles
export compiler
export cflags
export lflags
export libs
export includes

.PHONY: parties debug grid xdmf dry dryDebug clean cleanall

# Make 3-D code
parties: cflags += $(optmz_flags)
parties:
	$(MAKE) -C src

# Make 3-D code with debugging flags
debug: cflags += $(debug_flags) -DDEBUG
debug:
	$(MAKE) -C src

# Make grid generator
#     - generates grid for nonunifrom mesh
grid: cflags += $(optmz_flags)
grid:
	$(MAKE) -C grid_generator

# Make xdmf writer
#     - generates reader for HDF5 files
xdmf: cflags += $(optmz_flags)
xdmf:
	$(MAKE) -C xdmf_writer

# Make dry particle code
#     - simulates particles without fluid
dry: cflags += $(optmz_flags) -DDRY_PARTICLES
dry:
	$(MAKE) -C particles_dry

# Dry particle code with debugging
dryDebug: cflags += $(debug_flags) -DDRY_PARTICLES -DDEBUG
dryDebug:
	$(MAKE) -C particles_dry

# Make post-processing
#     - allows post-processing of HDF5 files
post: cflags += $(optmz_flags) -DPOST_PROCESS
post:
	$(MAKE) -C post_processing

# Make all targets
all: cart3d grid xdmf

clean:
	-rm -rf obj
	-rm -f parties
	-rm -f parties_dry

cleanall: clean
	-rm -f gridGenerator
	-rm -f xdmfWriter
	-rm -f postProcess

cleandata:
	-rm Data_*
	-rm Resume*
	-rm Reader_*.xmf
	-rm *.dat
	-rm *.h5
	-rm *.log
	-rm *.out
	-rm *.err
