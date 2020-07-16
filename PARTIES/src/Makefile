
# Add dependencies here
#-----------------------------------------------------------

# Files that all source files depend on -- These should use
# absolute pathname, or their location should be added to
# VPATH (below)
export DEPS := Boundary.h

# Add new subdirectories to this file
#-----------------------------------------------------------

# Defines SRC_* variables
include Makefile_src.def
SUB_DIRS := $(SRC_SUB_DIRS)
INCLUDES := $(SRC_INCLUDES)
OBJ_DIR  := $(SRC_OBJ_DIR)
OBJECTS  := $(SRC_OBJECTS)
ALL_OBJECTS := $(SRC_ALL_OBJECTS)

#-----------------------------------------------------------

# Include source code include files
includes += $(addprefix -I, $(INCLUDES))

# Location for object directory relative to here
export OBJ_ROOT_DIR := $(SRC_OBJ_ROOT_DIR)

# Make will use this path for searching for files
export VPATH := $(INCLUDES)

# Don't update just because timestamp on directories change
.PHONY: $(SUB_DIRS)

# Make objects in current directory and run makefiles in subdirectories
cart3d: $(OBJECTS) $(SUB_DIRS)
	$(compiler)  -o ../parties $(ALL_OBJECTS) $(lflags) $(libs) 

# Objects depend on existence of containing directory
$(OBJECTS): | $(OBJ_DIR)

# Create object directory, if it doesn't exist
$(OBJ_DIR):
	mkdir $(OBJ_DIR)

# Execute makefiles in subdirectories
$(SUB_DIRS):
	$(MAKE) -C $@ sources

# Compile object files based on local source files
$(OBJ_DIR)/%.o : %.c $(DEPS)
	$(compiler) $(cflags) $(includes)  -c $< -o $@
