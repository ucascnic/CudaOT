# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.26

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
CMAKE_COMMAND = /home/spack/spack/opt/spack/linux-debian12-zen2/gcc-11.3.0/cmake-3.26.3-gytjv5owjrxi7wdi3sjxvvcqf3nn26js/bin/cmake

# The command to remove a file.
RM = /home/spack/spack/opt/spack/linux-debian12-zen2/gcc-11.3.0/cmake-3.26.3-gytjv5owjrxi7wdi3sjxvvcqf3nn26js/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/chenyidong/CudaOT

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/chenyidong/CudaOT

# Include any dependencies generated for this target.
include ShortCutSolver/CMakeFiles/ShortCutSolver.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include ShortCutSolver/CMakeFiles/ShortCutSolver.dir/compiler_depend.make

# Include the progress variables for this target.
include ShortCutSolver/CMakeFiles/ShortCutSolver.dir/progress.make

# Include the compile flags for this target's objects.
include ShortCutSolver/CMakeFiles/ShortCutSolver.dir/flags.make

ShortCutSolver/CMakeFiles/ShortCutSolver.dir/Interfaces.cpp.o: ShortCutSolver/CMakeFiles/ShortCutSolver.dir/flags.make
ShortCutSolver/CMakeFiles/ShortCutSolver.dir/Interfaces.cpp.o: ShortCutSolver/Interfaces.cpp
ShortCutSolver/CMakeFiles/ShortCutSolver.dir/Interfaces.cpp.o: ShortCutSolver/CMakeFiles/ShortCutSolver.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/chenyidong/CudaOT/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object ShortCutSolver/CMakeFiles/ShortCutSolver.dir/Interfaces.cpp.o"
	cd /home/chenyidong/CudaOT/ShortCutSolver && /home/spack/spack/opt/spack/linux-debian11-zen2/gcc-10.2.1/gcc-8.5.0-ffkw5xy5ijm5nzbwjxyjpzs5usxxj3ur/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT ShortCutSolver/CMakeFiles/ShortCutSolver.dir/Interfaces.cpp.o -MF CMakeFiles/ShortCutSolver.dir/Interfaces.cpp.o.d -o CMakeFiles/ShortCutSolver.dir/Interfaces.cpp.o -c /home/chenyidong/CudaOT/ShortCutSolver/Interfaces.cpp

ShortCutSolver/CMakeFiles/ShortCutSolver.dir/Interfaces.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ShortCutSolver.dir/Interfaces.cpp.i"
	cd /home/chenyidong/CudaOT/ShortCutSolver && /home/spack/spack/opt/spack/linux-debian11-zen2/gcc-10.2.1/gcc-8.5.0-ffkw5xy5ijm5nzbwjxyjpzs5usxxj3ur/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/chenyidong/CudaOT/ShortCutSolver/Interfaces.cpp > CMakeFiles/ShortCutSolver.dir/Interfaces.cpp.i

ShortCutSolver/CMakeFiles/ShortCutSolver.dir/Interfaces.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ShortCutSolver.dir/Interfaces.cpp.s"
	cd /home/chenyidong/CudaOT/ShortCutSolver && /home/spack/spack/opt/spack/linux-debian11-zen2/gcc-10.2.1/gcc-8.5.0-ffkw5xy5ijm5nzbwjxyjpzs5usxxj3ur/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/chenyidong/CudaOT/ShortCutSolver/Interfaces.cpp -o CMakeFiles/ShortCutSolver.dir/Interfaces.cpp.s

ShortCutSolver/CMakeFiles/ShortCutSolver.dir/MultiScaleSolver.cpp.o: ShortCutSolver/CMakeFiles/ShortCutSolver.dir/flags.make
ShortCutSolver/CMakeFiles/ShortCutSolver.dir/MultiScaleSolver.cpp.o: ShortCutSolver/MultiScaleSolver.cpp
ShortCutSolver/CMakeFiles/ShortCutSolver.dir/MultiScaleSolver.cpp.o: ShortCutSolver/CMakeFiles/ShortCutSolver.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/chenyidong/CudaOT/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object ShortCutSolver/CMakeFiles/ShortCutSolver.dir/MultiScaleSolver.cpp.o"
	cd /home/chenyidong/CudaOT/ShortCutSolver && /home/spack/spack/opt/spack/linux-debian11-zen2/gcc-10.2.1/gcc-8.5.0-ffkw5xy5ijm5nzbwjxyjpzs5usxxj3ur/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT ShortCutSolver/CMakeFiles/ShortCutSolver.dir/MultiScaleSolver.cpp.o -MF CMakeFiles/ShortCutSolver.dir/MultiScaleSolver.cpp.o.d -o CMakeFiles/ShortCutSolver.dir/MultiScaleSolver.cpp.o -c /home/chenyidong/CudaOT/ShortCutSolver/MultiScaleSolver.cpp

ShortCutSolver/CMakeFiles/ShortCutSolver.dir/MultiScaleSolver.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ShortCutSolver.dir/MultiScaleSolver.cpp.i"
	cd /home/chenyidong/CudaOT/ShortCutSolver && /home/spack/spack/opt/spack/linux-debian11-zen2/gcc-10.2.1/gcc-8.5.0-ffkw5xy5ijm5nzbwjxyjpzs5usxxj3ur/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/chenyidong/CudaOT/ShortCutSolver/MultiScaleSolver.cpp > CMakeFiles/ShortCutSolver.dir/MultiScaleSolver.cpp.i

ShortCutSolver/CMakeFiles/ShortCutSolver.dir/MultiScaleSolver.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ShortCutSolver.dir/MultiScaleSolver.cpp.s"
	cd /home/chenyidong/CudaOT/ShortCutSolver && /home/spack/spack/opt/spack/linux-debian11-zen2/gcc-10.2.1/gcc-8.5.0-ffkw5xy5ijm5nzbwjxyjpzs5usxxj3ur/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/chenyidong/CudaOT/ShortCutSolver/MultiScaleSolver.cpp -o CMakeFiles/ShortCutSolver.dir/MultiScaleSolver.cpp.s

ShortCutSolver/CMakeFiles/ShortCutSolver.dir/TShortCutSolver.cpp.o: ShortCutSolver/CMakeFiles/ShortCutSolver.dir/flags.make
ShortCutSolver/CMakeFiles/ShortCutSolver.dir/TShortCutSolver.cpp.o: ShortCutSolver/TShortCutSolver.cpp
ShortCutSolver/CMakeFiles/ShortCutSolver.dir/TShortCutSolver.cpp.o: ShortCutSolver/CMakeFiles/ShortCutSolver.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/chenyidong/CudaOT/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object ShortCutSolver/CMakeFiles/ShortCutSolver.dir/TShortCutSolver.cpp.o"
	cd /home/chenyidong/CudaOT/ShortCutSolver && /home/spack/spack/opt/spack/linux-debian11-zen2/gcc-10.2.1/gcc-8.5.0-ffkw5xy5ijm5nzbwjxyjpzs5usxxj3ur/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT ShortCutSolver/CMakeFiles/ShortCutSolver.dir/TShortCutSolver.cpp.o -MF CMakeFiles/ShortCutSolver.dir/TShortCutSolver.cpp.o.d -o CMakeFiles/ShortCutSolver.dir/TShortCutSolver.cpp.o -c /home/chenyidong/CudaOT/ShortCutSolver/TShortCutSolver.cpp

ShortCutSolver/CMakeFiles/ShortCutSolver.dir/TShortCutSolver.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ShortCutSolver.dir/TShortCutSolver.cpp.i"
	cd /home/chenyidong/CudaOT/ShortCutSolver && /home/spack/spack/opt/spack/linux-debian11-zen2/gcc-10.2.1/gcc-8.5.0-ffkw5xy5ijm5nzbwjxyjpzs5usxxj3ur/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/chenyidong/CudaOT/ShortCutSolver/TShortCutSolver.cpp > CMakeFiles/ShortCutSolver.dir/TShortCutSolver.cpp.i

ShortCutSolver/CMakeFiles/ShortCutSolver.dir/TShortCutSolver.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ShortCutSolver.dir/TShortCutSolver.cpp.s"
	cd /home/chenyidong/CudaOT/ShortCutSolver && /home/spack/spack/opt/spack/linux-debian11-zen2/gcc-10.2.1/gcc-8.5.0-ffkw5xy5ijm5nzbwjxyjpzs5usxxj3ur/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/chenyidong/CudaOT/ShortCutSolver/TShortCutSolver.cpp -o CMakeFiles/ShortCutSolver.dir/TShortCutSolver.cpp.s

ShortCutSolver/CMakeFiles/ShortCutSolver.dir/TShieldGenerator.cpp.o: ShortCutSolver/CMakeFiles/ShortCutSolver.dir/flags.make
ShortCutSolver/CMakeFiles/ShortCutSolver.dir/TShieldGenerator.cpp.o: ShortCutSolver/TShieldGenerator.cpp
ShortCutSolver/CMakeFiles/ShortCutSolver.dir/TShieldGenerator.cpp.o: ShortCutSolver/CMakeFiles/ShortCutSolver.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/chenyidong/CudaOT/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object ShortCutSolver/CMakeFiles/ShortCutSolver.dir/TShieldGenerator.cpp.o"
	cd /home/chenyidong/CudaOT/ShortCutSolver && /home/spack/spack/opt/spack/linux-debian11-zen2/gcc-10.2.1/gcc-8.5.0-ffkw5xy5ijm5nzbwjxyjpzs5usxxj3ur/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT ShortCutSolver/CMakeFiles/ShortCutSolver.dir/TShieldGenerator.cpp.o -MF CMakeFiles/ShortCutSolver.dir/TShieldGenerator.cpp.o.d -o CMakeFiles/ShortCutSolver.dir/TShieldGenerator.cpp.o -c /home/chenyidong/CudaOT/ShortCutSolver/TShieldGenerator.cpp

ShortCutSolver/CMakeFiles/ShortCutSolver.dir/TShieldGenerator.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ShortCutSolver.dir/TShieldGenerator.cpp.i"
	cd /home/chenyidong/CudaOT/ShortCutSolver && /home/spack/spack/opt/spack/linux-debian11-zen2/gcc-10.2.1/gcc-8.5.0-ffkw5xy5ijm5nzbwjxyjpzs5usxxj3ur/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/chenyidong/CudaOT/ShortCutSolver/TShieldGenerator.cpp > CMakeFiles/ShortCutSolver.dir/TShieldGenerator.cpp.i

ShortCutSolver/CMakeFiles/ShortCutSolver.dir/TShieldGenerator.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ShortCutSolver.dir/TShieldGenerator.cpp.s"
	cd /home/chenyidong/CudaOT/ShortCutSolver && /home/spack/spack/opt/spack/linux-debian11-zen2/gcc-10.2.1/gcc-8.5.0-ffkw5xy5ijm5nzbwjxyjpzs5usxxj3ur/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/chenyidong/CudaOT/ShortCutSolver/TShieldGenerator.cpp -o CMakeFiles/ShortCutSolver.dir/TShieldGenerator.cpp.s

ShortCutSolver/CMakeFiles/ShortCutSolver.dir/TShieldGenerator_Models.cpp.o: ShortCutSolver/CMakeFiles/ShortCutSolver.dir/flags.make
ShortCutSolver/CMakeFiles/ShortCutSolver.dir/TShieldGenerator_Models.cpp.o: ShortCutSolver/TShieldGenerator_Models.cpp
ShortCutSolver/CMakeFiles/ShortCutSolver.dir/TShieldGenerator_Models.cpp.o: ShortCutSolver/CMakeFiles/ShortCutSolver.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/chenyidong/CudaOT/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object ShortCutSolver/CMakeFiles/ShortCutSolver.dir/TShieldGenerator_Models.cpp.o"
	cd /home/chenyidong/CudaOT/ShortCutSolver && /home/spack/spack/opt/spack/linux-debian11-zen2/gcc-10.2.1/gcc-8.5.0-ffkw5xy5ijm5nzbwjxyjpzs5usxxj3ur/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT ShortCutSolver/CMakeFiles/ShortCutSolver.dir/TShieldGenerator_Models.cpp.o -MF CMakeFiles/ShortCutSolver.dir/TShieldGenerator_Models.cpp.o.d -o CMakeFiles/ShortCutSolver.dir/TShieldGenerator_Models.cpp.o -c /home/chenyidong/CudaOT/ShortCutSolver/TShieldGenerator_Models.cpp

ShortCutSolver/CMakeFiles/ShortCutSolver.dir/TShieldGenerator_Models.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ShortCutSolver.dir/TShieldGenerator_Models.cpp.i"
	cd /home/chenyidong/CudaOT/ShortCutSolver && /home/spack/spack/opt/spack/linux-debian11-zen2/gcc-10.2.1/gcc-8.5.0-ffkw5xy5ijm5nzbwjxyjpzs5usxxj3ur/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/chenyidong/CudaOT/ShortCutSolver/TShieldGenerator_Models.cpp > CMakeFiles/ShortCutSolver.dir/TShieldGenerator_Models.cpp.i

ShortCutSolver/CMakeFiles/ShortCutSolver.dir/TShieldGenerator_Models.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ShortCutSolver.dir/TShieldGenerator_Models.cpp.s"
	cd /home/chenyidong/CudaOT/ShortCutSolver && /home/spack/spack/opt/spack/linux-debian11-zen2/gcc-10.2.1/gcc-8.5.0-ffkw5xy5ijm5nzbwjxyjpzs5usxxj3ur/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/chenyidong/CudaOT/ShortCutSolver/TShieldGenerator_Models.cpp -o CMakeFiles/ShortCutSolver.dir/TShieldGenerator_Models.cpp.s

# Object files for target ShortCutSolver
ShortCutSolver_OBJECTS = \
"CMakeFiles/ShortCutSolver.dir/Interfaces.cpp.o" \
"CMakeFiles/ShortCutSolver.dir/MultiScaleSolver.cpp.o" \
"CMakeFiles/ShortCutSolver.dir/TShortCutSolver.cpp.o" \
"CMakeFiles/ShortCutSolver.dir/TShieldGenerator.cpp.o" \
"CMakeFiles/ShortCutSolver.dir/TShieldGenerator_Models.cpp.o"

# External object files for target ShortCutSolver
ShortCutSolver_EXTERNAL_OBJECTS =

ShortCutSolver/libShortCutSolver.a: ShortCutSolver/CMakeFiles/ShortCutSolver.dir/Interfaces.cpp.o
ShortCutSolver/libShortCutSolver.a: ShortCutSolver/CMakeFiles/ShortCutSolver.dir/MultiScaleSolver.cpp.o
ShortCutSolver/libShortCutSolver.a: ShortCutSolver/CMakeFiles/ShortCutSolver.dir/TShortCutSolver.cpp.o
ShortCutSolver/libShortCutSolver.a: ShortCutSolver/CMakeFiles/ShortCutSolver.dir/TShieldGenerator.cpp.o
ShortCutSolver/libShortCutSolver.a: ShortCutSolver/CMakeFiles/ShortCutSolver.dir/TShieldGenerator_Models.cpp.o
ShortCutSolver/libShortCutSolver.a: ShortCutSolver/CMakeFiles/ShortCutSolver.dir/build.make
ShortCutSolver/libShortCutSolver.a: ShortCutSolver/CMakeFiles/ShortCutSolver.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/chenyidong/CudaOT/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Linking CXX static library libShortCutSolver.a"
	cd /home/chenyidong/CudaOT/ShortCutSolver && $(CMAKE_COMMAND) -P CMakeFiles/ShortCutSolver.dir/cmake_clean_target.cmake
	cd /home/chenyidong/CudaOT/ShortCutSolver && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/ShortCutSolver.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
ShortCutSolver/CMakeFiles/ShortCutSolver.dir/build: ShortCutSolver/libShortCutSolver.a
.PHONY : ShortCutSolver/CMakeFiles/ShortCutSolver.dir/build

ShortCutSolver/CMakeFiles/ShortCutSolver.dir/clean:
	cd /home/chenyidong/CudaOT/ShortCutSolver && $(CMAKE_COMMAND) -P CMakeFiles/ShortCutSolver.dir/cmake_clean.cmake
.PHONY : ShortCutSolver/CMakeFiles/ShortCutSolver.dir/clean

ShortCutSolver/CMakeFiles/ShortCutSolver.dir/depend:
	cd /home/chenyidong/CudaOT && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/chenyidong/CudaOT /home/chenyidong/CudaOT/ShortCutSolver /home/chenyidong/CudaOT /home/chenyidong/CudaOT/ShortCutSolver /home/chenyidong/CudaOT/ShortCutSolver/CMakeFiles/ShortCutSolver.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : ShortCutSolver/CMakeFiles/ShortCutSolver.dir/depend
