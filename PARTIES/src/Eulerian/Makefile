
# Local directory
#-----------------------------------------------------------
LOCAL_DIR := $(notdir $(CURDIR))
OBJ_DIR   := ../$(OBJ_ROOT_DIR)/$(LOCAL_DIR)

# Add dependencies here
#
# E.g. Pressure.c depends on 'lsolver/psolve_fft.c':
#
# $(OBJ_DIR)/Pressure.o : lsolver/psolve_fft.c
#-----------------------------------------------------------

$(OBJ_DIR)/Velocity.o : lsolver/msolve_direct.c lsolver/msolve_cg.c lsolver/msolve_bicg.c lsolver/direct_solver.c

$(OBJ_DIR)/Pressure.o : lsolver/psolve_fft.c

$(OBJ_DIR)/Conc.o : lsolver/csolve_direct.c

# Do not edit below here
#***********************************************************

# Source and object files associated with this directory
SOURCES = $(wildcard *.c)
OBJECTS = $(patsubst %.c, $(OBJ_DIR)/%.o, $(SOURCES))

# Make objects in this directory
sources: $(OBJECTS) $(OBJ_DIR)

# Objects depend on existence of containing directory
$(OBJECTS): | $(OBJ_DIR)

# Create object directory, if it doesn't exist
$(OBJ_DIR):
	mkdir $(OBJ_DIR)

# Compile object files based on local source files
$(OBJ_DIR)/%.o : %.c $(DEPS)
	$(compiler) $(cflags) $(includes) -c $< -o $@
