# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.24

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
CMAKE_COMMAND = /usr/local/bin/cmake

# The command to remove a file.
RM = /usr/local/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/zhyu/projects/cna_clone/software/prep

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/zhyu/projects/cna_clone/software/prep

# Include any dependencies generated for this target.
include lib/CMakeFiles/split.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include lib/CMakeFiles/split.dir/compiler_depend.make

# Include the progress variables for this target.
include lib/CMakeFiles/split.dir/progress.make

# Include the compile flags for this target's objects.
include lib/CMakeFiles/split.dir/flags.make

lib/CMakeFiles/split.dir/split/split.cpp.o: lib/CMakeFiles/split.dir/flags.make
lib/CMakeFiles/split.dir/split/split.cpp.o: lib/split/split.cpp
lib/CMakeFiles/split.dir/split/split.cpp.o: lib/CMakeFiles/split.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/zhyu/projects/cna_clone/software/prep/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object lib/CMakeFiles/split.dir/split/split.cpp.o"
	cd /home/zhyu/projects/cna_clone/software/prep/lib && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT lib/CMakeFiles/split.dir/split/split.cpp.o -MF CMakeFiles/split.dir/split/split.cpp.o.d -o CMakeFiles/split.dir/split/split.cpp.o -c /home/zhyu/projects/cna_clone/software/prep/lib/split/split.cpp

lib/CMakeFiles/split.dir/split/split.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/split.dir/split/split.cpp.i"
	cd /home/zhyu/projects/cna_clone/software/prep/lib && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/zhyu/projects/cna_clone/software/prep/lib/split/split.cpp > CMakeFiles/split.dir/split/split.cpp.i

lib/CMakeFiles/split.dir/split/split.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/split.dir/split/split.cpp.s"
	cd /home/zhyu/projects/cna_clone/software/prep/lib && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/zhyu/projects/cna_clone/software/prep/lib/split/split.cpp -o CMakeFiles/split.dir/split/split.cpp.s

# Object files for target split
split_OBJECTS = \
"CMakeFiles/split.dir/split/split.cpp.o"

# External object files for target split
split_EXTERNAL_OBJECTS =

lib/libsplit.a: lib/CMakeFiles/split.dir/split/split.cpp.o
lib/libsplit.a: lib/CMakeFiles/split.dir/build.make
lib/libsplit.a: lib/CMakeFiles/split.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/zhyu/projects/cna_clone/software/prep/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX static library libsplit.a"
	cd /home/zhyu/projects/cna_clone/software/prep/lib && $(CMAKE_COMMAND) -P CMakeFiles/split.dir/cmake_clean_target.cmake
	cd /home/zhyu/projects/cna_clone/software/prep/lib && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/split.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
lib/CMakeFiles/split.dir/build: lib/libsplit.a
.PHONY : lib/CMakeFiles/split.dir/build

lib/CMakeFiles/split.dir/clean:
	cd /home/zhyu/projects/cna_clone/software/prep/lib && $(CMAKE_COMMAND) -P CMakeFiles/split.dir/cmake_clean.cmake
.PHONY : lib/CMakeFiles/split.dir/clean

lib/CMakeFiles/split.dir/depend:
	cd /home/zhyu/projects/cna_clone/software/prep && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/zhyu/projects/cna_clone/software/prep /home/zhyu/projects/cna_clone/software/prep/lib /home/zhyu/projects/cna_clone/software/prep /home/zhyu/projects/cna_clone/software/prep/lib /home/zhyu/projects/cna_clone/software/prep/lib/CMakeFiles/split.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : lib/CMakeFiles/split.dir/depend
