# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.20

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /snap/clion/169/bin/cmake/linux/bin/cmake

# The command to remove a file.
RM = /snap/clion/169/bin/cmake/linux/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/sergio/repos/gravity_simulator_omp

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/sergio/repos/gravity_simulator_omp/cmake-build-release

# Include any dependencies generated for this target.
include CMakeFiles/paos.dir/depend.make
# Include the progress variables for this target.
include CMakeFiles/paos.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/paos.dir/flags.make

CMakeFiles/paos.dir/paos.cpp.o: CMakeFiles/paos.dir/flags.make
CMakeFiles/paos.dir/paos.cpp.o: ../paos.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/sergio/repos/gravity_simulator_omp/cmake-build-release/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/paos.dir/paos.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/paos.dir/paos.cpp.o -c /home/sergio/repos/gravity_simulator_omp/paos.cpp

CMakeFiles/paos.dir/paos.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/paos.dir/paos.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/sergio/repos/gravity_simulator_omp/paos.cpp > CMakeFiles/paos.dir/paos.cpp.i

CMakeFiles/paos.dir/paos.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/paos.dir/paos.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/sergio/repos/gravity_simulator_omp/paos.cpp -o CMakeFiles/paos.dir/paos.cpp.s

# Object files for target paos
paos_OBJECTS = \
"CMakeFiles/paos.dir/paos.cpp.o"

# External object files for target paos
paos_EXTERNAL_OBJECTS =

paos: CMakeFiles/paos.dir/paos.cpp.o
paos: CMakeFiles/paos.dir/build.make
paos: /usr/lib/gcc/x86_64-linux-gnu/9/libgomp.so
paos: /usr/lib/x86_64-linux-gnu/libpthread.so
paos: CMakeFiles/paos.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/sergio/repos/gravity_simulator_omp/cmake-build-release/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable paos"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/paos.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/paos.dir/build: paos
.PHONY : CMakeFiles/paos.dir/build

CMakeFiles/paos.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/paos.dir/cmake_clean.cmake
.PHONY : CMakeFiles/paos.dir/clean

CMakeFiles/paos.dir/depend:
	cd /home/sergio/repos/gravity_simulator_omp/cmake-build-release && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/sergio/repos/gravity_simulator_omp /home/sergio/repos/gravity_simulator_omp /home/sergio/repos/gravity_simulator_omp/cmake-build-release /home/sergio/repos/gravity_simulator_omp/cmake-build-release /home/sergio/repos/gravity_simulator_omp/cmake-build-release/CMakeFiles/paos.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/paos.dir/depend

