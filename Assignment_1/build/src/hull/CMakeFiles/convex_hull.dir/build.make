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
include src/hull/CMakeFiles/convex_hull.dir/depend.make

# Include the progress variables for this target.
include src/hull/CMakeFiles/convex_hull.dir/progress.make

# Include the compile flags for this target's objects.
include src/hull/CMakeFiles/convex_hull.dir/flags.make

src/hull/CMakeFiles/convex_hull.dir/main.cpp.o: src/hull/CMakeFiles/convex_hull.dir/flags.make
src/hull/CMakeFiles/convex_hull.dir/main.cpp.o: ../src/hull/main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/zdrm/Desktop/Fall19/CG/ucg/Assignment_1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/hull/CMakeFiles/convex_hull.dir/main.cpp.o"
	cd /home/zdrm/Desktop/Fall19/CG/ucg/Assignment_1/build/src/hull && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/convex_hull.dir/main.cpp.o -c /home/zdrm/Desktop/Fall19/CG/ucg/Assignment_1/src/hull/main.cpp

src/hull/CMakeFiles/convex_hull.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/convex_hull.dir/main.cpp.i"
	cd /home/zdrm/Desktop/Fall19/CG/ucg/Assignment_1/build/src/hull && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/zdrm/Desktop/Fall19/CG/ucg/Assignment_1/src/hull/main.cpp > CMakeFiles/convex_hull.dir/main.cpp.i

src/hull/CMakeFiles/convex_hull.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/convex_hull.dir/main.cpp.s"
	cd /home/zdrm/Desktop/Fall19/CG/ucg/Assignment_1/build/src/hull && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/zdrm/Desktop/Fall19/CG/ucg/Assignment_1/src/hull/main.cpp -o CMakeFiles/convex_hull.dir/main.cpp.s

src/hull/CMakeFiles/convex_hull.dir/main.cpp.o.requires:

.PHONY : src/hull/CMakeFiles/convex_hull.dir/main.cpp.o.requires

src/hull/CMakeFiles/convex_hull.dir/main.cpp.o.provides: src/hull/CMakeFiles/convex_hull.dir/main.cpp.o.requires
	$(MAKE) -f src/hull/CMakeFiles/convex_hull.dir/build.make src/hull/CMakeFiles/convex_hull.dir/main.cpp.o.provides.build
.PHONY : src/hull/CMakeFiles/convex_hull.dir/main.cpp.o.provides

src/hull/CMakeFiles/convex_hull.dir/main.cpp.o.provides.build: src/hull/CMakeFiles/convex_hull.dir/main.cpp.o


# Object files for target convex_hull
convex_hull_OBJECTS = \
"CMakeFiles/convex_hull.dir/main.cpp.o"

# External object files for target convex_hull
convex_hull_EXTERNAL_OBJECTS =

convex_hull: src/hull/CMakeFiles/convex_hull.dir/main.cpp.o
convex_hull: src/hull/CMakeFiles/convex_hull.dir/build.make
convex_hull: src/hull/CMakeFiles/convex_hull.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/zdrm/Desktop/Fall19/CG/ucg/Assignment_1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../../convex_hull"
	cd /home/zdrm/Desktop/Fall19/CG/ucg/Assignment_1/build/src/hull && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/convex_hull.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/hull/CMakeFiles/convex_hull.dir/build: convex_hull

.PHONY : src/hull/CMakeFiles/convex_hull.dir/build

src/hull/CMakeFiles/convex_hull.dir/requires: src/hull/CMakeFiles/convex_hull.dir/main.cpp.o.requires

.PHONY : src/hull/CMakeFiles/convex_hull.dir/requires

src/hull/CMakeFiles/convex_hull.dir/clean:
	cd /home/zdrm/Desktop/Fall19/CG/ucg/Assignment_1/build/src/hull && $(CMAKE_COMMAND) -P CMakeFiles/convex_hull.dir/cmake_clean.cmake
.PHONY : src/hull/CMakeFiles/convex_hull.dir/clean

src/hull/CMakeFiles/convex_hull.dir/depend:
	cd /home/zdrm/Desktop/Fall19/CG/ucg/Assignment_1/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/zdrm/Desktop/Fall19/CG/ucg/Assignment_1 /home/zdrm/Desktop/Fall19/CG/ucg/Assignment_1/src/hull /home/zdrm/Desktop/Fall19/CG/ucg/Assignment_1/build /home/zdrm/Desktop/Fall19/CG/ucg/Assignment_1/build/src/hull /home/zdrm/Desktop/Fall19/CG/ucg/Assignment_1/build/src/hull/CMakeFiles/convex_hull.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/hull/CMakeFiles/convex_hull.dir/depend

