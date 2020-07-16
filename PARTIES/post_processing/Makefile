
# Get variables SRC_INCLUDES and SRC_ALL_OBJECTS
include ../src/Makefile_src.def

# Filter out object files that contain main() functions
SRC_ALL_OBJECTS := $(filter-out ../obj/Cart3d.o, $(SRC_ALL_OBJECTS))

# Local source files
SOURCES := $(wildcard *.c)

# Make sure to run cart3d recipe
.PHONY: cart3d

# Compile local files after cart3d objects have been made
post: includes += $(addprefix -I, $(SRC_INCLUDES))
post: cart3d
	$(compiler) $(cflags) $(includes) $(SOURCES) -o ../postProcess $(SRC_ALL_OBJECTS) $(lflags) $(libs)

# Make all object files for cart3d
cart3d:
	$(MAKE) -C ../src cart3d
