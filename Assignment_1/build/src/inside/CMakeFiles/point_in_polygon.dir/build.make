# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/zdrm/Desktop/Fall19/CG/ucg/Assignment_1

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/zdrm/Desktop/Fall19/CG/ucg/Assignment_1/build

# Include any dependencies generated for this target.
include src/inside/CMakeFiles/point_in_polygon.dir/depend.make

# Include the progress variables for this target.
include src/inside/CMakeFiles/point_in_polygon.dir/progress.make

# Include the compile flags for this target's objects.
include src/inside/CMakeFiles/point_in_polygon.dir/flags.make

src/inside/CMakeFiles/point_in_polygon.dir/main.cpp.o: src/inside/CMakeFiles/point_in_polygon.dir/flags.make
src/inside/CMakeFiles/point_in_polygon.dir/main.cpp.o: ../src/inside/main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/zdrm/Desktop/Fall19/CG/ucg/Assignment_1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/inside/CMakeFiles/point_in_polygon.dir/main.cpp.o"
	cd /home/zdrm/Desktop/Fall19/CG/ucg/Assignment_1/build/src/inside && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/point_in_polygon.dir/main.cpp.o -c /home/zdrm/Desktop/Fall19/CG/ucg/Assignment_1/src/inside/main.cpp

src/inside/CMakeFiles/point_in_polygon.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/point_in_polygon.dir/main.cpp.i"
	cd /home/zdrm/Desktop/Fall19/CG/ucg/Assignment_1/build/src/inside && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/zdrm/Desktop/Fall19/CG/ucg/Assignment_1/src/inside/main.cpp > CMakeFiles/point_in_polygon.dir/main.cpp.i

src/inside/CMakeFiles/point_in_polygon.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/point_in_polygon.dir/main.cpp.s"
	cd /home/zdrm/Desktop/Fall19/CG/ucg/Assignment_1/build/src/inside && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/zdrm/Desktop/Fall19/CG/ucg/Assignment_1/src/inside/main.cpp -o CMakeFiles/point_in_polygon.dir/main.cpp.s

src/inside/CMakeFiles/point_in_polygon.dir/main.cpp.o.requires:

.PHONY : src/inside/CMakeFiles/point_in_polygon.dir/main.cpp.o.requires

src/inside/CMakeFiles/point_in_polygon.dir/main.cpp.o.provides: src/inside/CMakeFiles/point_in_polygon.dir/main.cpp.o.requires
	$(MAKE) -f src/inside/CMakeFiles/point_in_polygon.dir/build.make src/inside/CMakeFiles/point_in_polygon.dir/main.cpp.o.provides.build
.PHONY : src/inside/CMakeFiles/point_in_polygon.dir/main.cpp.o.provides

src/inside/CMakeFiles/point_in_polygon.dir/main.cpp.o.provides.build: src/inside/CMakeFiles/point_in_polygon.dir/main.cpp.o


# Object files for target point_in_polygon
point_in_polygon_OBJECTS = \
"CMakeFiles/point_in_polygon.dir/main.cpp.o"

# External object files for target point_in_polygon
point_in_polygon_EXTERNAL_OBJECTS =

point_in_polygon: src/inside/CMakeFiles/point_in_polygon.dir/main.cpp.o
point_in_polygon: src/inside/CMakeFiles/point_in_polygon.dir/build.make
point_in_polygon: src/inside/CMakeFiles/point_in_polygon.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/zdrm/Desktop/Fall19/CG/ucg/Assignment_1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../../point_in_polygon"
	cd /home/zdrm/Desktop/Fall19/CG/ucg/Assignment_1/build/src/inside && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/point_in_polygon.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/inside/CMakeFiles/point_in_polygon.dir/build: point_in_polygon

.PHONY : src/inside/CMakeFiles/point_in_polygon.dir/build

src/inside/CMakeFiles/point_in_polygon.dir/requires: src/inside/CMakeFiles/point_in_polygon.dir/main.cpp.o.requires

.PHONY : src/inside/CMakeFiles/point_in_polygon.dir/requires

src/inside/CMakeFiles/point_in_polygon.dir/clean:
	cd /home/zdrm/Desktop/Fall19/CG/ucg/Assignment_1/build/src/inside && $(CMAKE_COMMAND) -P CMakeFiles/point_in_polygon.dir/cmake_clean.cmake
.PHONY : src/inside/CMakeFiles/point_in_polygon.dir/clean

src/inside/CMakeFiles/point_in_polygon.dir/depend:
	cd /home/zdrm/Desktop/Fall19/CG/ucg/Assignment_1/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/zdrm/Desktop/Fall19/CG/ucg/Assignment_1 /home/zdrm/Desktop/Fall19/CG/ucg/Assignment_1/src/inside /home/zdrm/Desktop/Fall19/CG/ucg/Assignment_1/build /home/zdrm/Desktop/Fall19/CG/ucg/Assignment_1/build/src/inside /home/zdrm/Desktop/Fall19/CG/ucg/Assignment_1/build/src/inside/CMakeFiles/point_in_polygon.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/inside/CMakeFiles/point_in_polygon.dir/depend

