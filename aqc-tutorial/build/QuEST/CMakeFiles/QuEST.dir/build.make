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
CMAKE_COMMAND = /usr/local/Cellar/cmake/3.26.0/bin/cmake

# The command to remove a file.
RM = /usr/local/Cellar/cmake/3.26.0/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/cathan/PhD/事务-课题组/Tutorial-website/2-sat

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/cathan/PhD/事务-课题组/Tutorial-website/2-sat/build

# Include any dependencies generated for this target.
include QuEST/CMakeFiles/QuEST.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include QuEST/CMakeFiles/QuEST.dir/compiler_depend.make

# Include the progress variables for this target.
include QuEST/CMakeFiles/QuEST.dir/progress.make

# Include the compile flags for this target's objects.
include QuEST/CMakeFiles/QuEST.dir/flags.make

QuEST/CMakeFiles/QuEST.dir/src/QuEST.c.o: QuEST/CMakeFiles/QuEST.dir/flags.make
QuEST/CMakeFiles/QuEST.dir/src/QuEST.c.o: /Users/cathan/PhD/事务-课题组/Tutorial-website/2-sat/QuEST/src/QuEST.c
QuEST/CMakeFiles/QuEST.dir/src/QuEST.c.o: QuEST/CMakeFiles/QuEST.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/cathan/PhD/事务-课题组/Tutorial-website/2-sat/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object QuEST/CMakeFiles/QuEST.dir/src/QuEST.c.o"
	cd /Users/cathan/PhD/事务-课题组/Tutorial-website/2-sat/build/QuEST && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT QuEST/CMakeFiles/QuEST.dir/src/QuEST.c.o -MF CMakeFiles/QuEST.dir/src/QuEST.c.o.d -o CMakeFiles/QuEST.dir/src/QuEST.c.o -c /Users/cathan/PhD/事务-课题组/Tutorial-website/2-sat/QuEST/src/QuEST.c

QuEST/CMakeFiles/QuEST.dir/src/QuEST.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/QuEST.dir/src/QuEST.c.i"
	cd /Users/cathan/PhD/事务-课题组/Tutorial-website/2-sat/build/QuEST && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/cathan/PhD/事务-课题组/Tutorial-website/2-sat/QuEST/src/QuEST.c > CMakeFiles/QuEST.dir/src/QuEST.c.i

QuEST/CMakeFiles/QuEST.dir/src/QuEST.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/QuEST.dir/src/QuEST.c.s"
	cd /Users/cathan/PhD/事务-课题组/Tutorial-website/2-sat/build/QuEST && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/cathan/PhD/事务-课题组/Tutorial-website/2-sat/QuEST/src/QuEST.c -o CMakeFiles/QuEST.dir/src/QuEST.c.s

QuEST/CMakeFiles/QuEST.dir/src/QuEST_common.c.o: QuEST/CMakeFiles/QuEST.dir/flags.make
QuEST/CMakeFiles/QuEST.dir/src/QuEST_common.c.o: /Users/cathan/PhD/事务-课题组/Tutorial-website/2-sat/QuEST/src/QuEST_common.c
QuEST/CMakeFiles/QuEST.dir/src/QuEST_common.c.o: QuEST/CMakeFiles/QuEST.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/cathan/PhD/事务-课题组/Tutorial-website/2-sat/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building C object QuEST/CMakeFiles/QuEST.dir/src/QuEST_common.c.o"
	cd /Users/cathan/PhD/事务-课题组/Tutorial-website/2-sat/build/QuEST && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT QuEST/CMakeFiles/QuEST.dir/src/QuEST_common.c.o -MF CMakeFiles/QuEST.dir/src/QuEST_common.c.o.d -o CMakeFiles/QuEST.dir/src/QuEST_common.c.o -c /Users/cathan/PhD/事务-课题组/Tutorial-website/2-sat/QuEST/src/QuEST_common.c

QuEST/CMakeFiles/QuEST.dir/src/QuEST_common.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/QuEST.dir/src/QuEST_common.c.i"
	cd /Users/cathan/PhD/事务-课题组/Tutorial-website/2-sat/build/QuEST && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/cathan/PhD/事务-课题组/Tutorial-website/2-sat/QuEST/src/QuEST_common.c > CMakeFiles/QuEST.dir/src/QuEST_common.c.i

QuEST/CMakeFiles/QuEST.dir/src/QuEST_common.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/QuEST.dir/src/QuEST_common.c.s"
	cd /Users/cathan/PhD/事务-课题组/Tutorial-website/2-sat/build/QuEST && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/cathan/PhD/事务-课题组/Tutorial-website/2-sat/QuEST/src/QuEST_common.c -o CMakeFiles/QuEST.dir/src/QuEST_common.c.s

QuEST/CMakeFiles/QuEST.dir/src/QuEST_qasm.c.o: QuEST/CMakeFiles/QuEST.dir/flags.make
QuEST/CMakeFiles/QuEST.dir/src/QuEST_qasm.c.o: /Users/cathan/PhD/事务-课题组/Tutorial-website/2-sat/QuEST/src/QuEST_qasm.c
QuEST/CMakeFiles/QuEST.dir/src/QuEST_qasm.c.o: QuEST/CMakeFiles/QuEST.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/cathan/PhD/事务-课题组/Tutorial-website/2-sat/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building C object QuEST/CMakeFiles/QuEST.dir/src/QuEST_qasm.c.o"
	cd /Users/cathan/PhD/事务-课题组/Tutorial-website/2-sat/build/QuEST && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT QuEST/CMakeFiles/QuEST.dir/src/QuEST_qasm.c.o -MF CMakeFiles/QuEST.dir/src/QuEST_qasm.c.o.d -o CMakeFiles/QuEST.dir/src/QuEST_qasm.c.o -c /Users/cathan/PhD/事务-课题组/Tutorial-website/2-sat/QuEST/src/QuEST_qasm.c

QuEST/CMakeFiles/QuEST.dir/src/QuEST_qasm.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/QuEST.dir/src/QuEST_qasm.c.i"
	cd /Users/cathan/PhD/事务-课题组/Tutorial-website/2-sat/build/QuEST && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/cathan/PhD/事务-课题组/Tutorial-website/2-sat/QuEST/src/QuEST_qasm.c > CMakeFiles/QuEST.dir/src/QuEST_qasm.c.i

QuEST/CMakeFiles/QuEST.dir/src/QuEST_qasm.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/QuEST.dir/src/QuEST_qasm.c.s"
	cd /Users/cathan/PhD/事务-课题组/Tutorial-website/2-sat/build/QuEST && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/cathan/PhD/事务-课题组/Tutorial-website/2-sat/QuEST/src/QuEST_qasm.c -o CMakeFiles/QuEST.dir/src/QuEST_qasm.c.s

QuEST/CMakeFiles/QuEST.dir/src/QuEST_validation.c.o: QuEST/CMakeFiles/QuEST.dir/flags.make
QuEST/CMakeFiles/QuEST.dir/src/QuEST_validation.c.o: /Users/cathan/PhD/事务-课题组/Tutorial-website/2-sat/QuEST/src/QuEST_validation.c
QuEST/CMakeFiles/QuEST.dir/src/QuEST_validation.c.o: QuEST/CMakeFiles/QuEST.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/cathan/PhD/事务-课题组/Tutorial-website/2-sat/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building C object QuEST/CMakeFiles/QuEST.dir/src/QuEST_validation.c.o"
	cd /Users/cathan/PhD/事务-课题组/Tutorial-website/2-sat/build/QuEST && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT QuEST/CMakeFiles/QuEST.dir/src/QuEST_validation.c.o -MF CMakeFiles/QuEST.dir/src/QuEST_validation.c.o.d -o CMakeFiles/QuEST.dir/src/QuEST_validation.c.o -c /Users/cathan/PhD/事务-课题组/Tutorial-website/2-sat/QuEST/src/QuEST_validation.c

QuEST/CMakeFiles/QuEST.dir/src/QuEST_validation.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/QuEST.dir/src/QuEST_validation.c.i"
	cd /Users/cathan/PhD/事务-课题组/Tutorial-website/2-sat/build/QuEST && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/cathan/PhD/事务-课题组/Tutorial-website/2-sat/QuEST/src/QuEST_validation.c > CMakeFiles/QuEST.dir/src/QuEST_validation.c.i

QuEST/CMakeFiles/QuEST.dir/src/QuEST_validation.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/QuEST.dir/src/QuEST_validation.c.s"
	cd /Users/cathan/PhD/事务-课题组/Tutorial-website/2-sat/build/QuEST && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/cathan/PhD/事务-课题组/Tutorial-website/2-sat/QuEST/src/QuEST_validation.c -o CMakeFiles/QuEST.dir/src/QuEST_validation.c.s

QuEST/CMakeFiles/QuEST.dir/src/mt19937ar.c.o: QuEST/CMakeFiles/QuEST.dir/flags.make
QuEST/CMakeFiles/QuEST.dir/src/mt19937ar.c.o: /Users/cathan/PhD/事务-课题组/Tutorial-website/2-sat/QuEST/src/mt19937ar.c
QuEST/CMakeFiles/QuEST.dir/src/mt19937ar.c.o: QuEST/CMakeFiles/QuEST.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/cathan/PhD/事务-课题组/Tutorial-website/2-sat/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building C object QuEST/CMakeFiles/QuEST.dir/src/mt19937ar.c.o"
	cd /Users/cathan/PhD/事务-课题组/Tutorial-website/2-sat/build/QuEST && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT QuEST/CMakeFiles/QuEST.dir/src/mt19937ar.c.o -MF CMakeFiles/QuEST.dir/src/mt19937ar.c.o.d -o CMakeFiles/QuEST.dir/src/mt19937ar.c.o -c /Users/cathan/PhD/事务-课题组/Tutorial-website/2-sat/QuEST/src/mt19937ar.c

QuEST/CMakeFiles/QuEST.dir/src/mt19937ar.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/QuEST.dir/src/mt19937ar.c.i"
	cd /Users/cathan/PhD/事务-课题组/Tutorial-website/2-sat/build/QuEST && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/cathan/PhD/事务-课题组/Tutorial-website/2-sat/QuEST/src/mt19937ar.c > CMakeFiles/QuEST.dir/src/mt19937ar.c.i

QuEST/CMakeFiles/QuEST.dir/src/mt19937ar.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/QuEST.dir/src/mt19937ar.c.s"
	cd /Users/cathan/PhD/事务-课题组/Tutorial-website/2-sat/build/QuEST && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/cathan/PhD/事务-课题组/Tutorial-website/2-sat/QuEST/src/mt19937ar.c -o CMakeFiles/QuEST.dir/src/mt19937ar.c.s

QuEST/CMakeFiles/QuEST.dir/src/CPU/QuEST_cpu.c.o: QuEST/CMakeFiles/QuEST.dir/flags.make
QuEST/CMakeFiles/QuEST.dir/src/CPU/QuEST_cpu.c.o: /Users/cathan/PhD/事务-课题组/Tutorial-website/2-sat/QuEST/src/CPU/QuEST_cpu.c
QuEST/CMakeFiles/QuEST.dir/src/CPU/QuEST_cpu.c.o: QuEST/CMakeFiles/QuEST.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/cathan/PhD/事务-课题组/Tutorial-website/2-sat/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building C object QuEST/CMakeFiles/QuEST.dir/src/CPU/QuEST_cpu.c.o"
	cd /Users/cathan/PhD/事务-课题组/Tutorial-website/2-sat/build/QuEST && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT QuEST/CMakeFiles/QuEST.dir/src/CPU/QuEST_cpu.c.o -MF CMakeFiles/QuEST.dir/src/CPU/QuEST_cpu.c.o.d -o CMakeFiles/QuEST.dir/src/CPU/QuEST_cpu.c.o -c /Users/cathan/PhD/事务-课题组/Tutorial-website/2-sat/QuEST/src/CPU/QuEST_cpu.c

QuEST/CMakeFiles/QuEST.dir/src/CPU/QuEST_cpu.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/QuEST.dir/src/CPU/QuEST_cpu.c.i"
	cd /Users/cathan/PhD/事务-课题组/Tutorial-website/2-sat/build/QuEST && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/cathan/PhD/事务-课题组/Tutorial-website/2-sat/QuEST/src/CPU/QuEST_cpu.c > CMakeFiles/QuEST.dir/src/CPU/QuEST_cpu.c.i

QuEST/CMakeFiles/QuEST.dir/src/CPU/QuEST_cpu.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/QuEST.dir/src/CPU/QuEST_cpu.c.s"
	cd /Users/cathan/PhD/事务-课题组/Tutorial-website/2-sat/build/QuEST && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/cathan/PhD/事务-课题组/Tutorial-website/2-sat/QuEST/src/CPU/QuEST_cpu.c -o CMakeFiles/QuEST.dir/src/CPU/QuEST_cpu.c.s

QuEST/CMakeFiles/QuEST.dir/src/CPU/QuEST_cpu_local.c.o: QuEST/CMakeFiles/QuEST.dir/flags.make
QuEST/CMakeFiles/QuEST.dir/src/CPU/QuEST_cpu_local.c.o: /Users/cathan/PhD/事务-课题组/Tutorial-website/2-sat/QuEST/src/CPU/QuEST_cpu_local.c
QuEST/CMakeFiles/QuEST.dir/src/CPU/QuEST_cpu_local.c.o: QuEST/CMakeFiles/QuEST.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/cathan/PhD/事务-课题组/Tutorial-website/2-sat/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building C object QuEST/CMakeFiles/QuEST.dir/src/CPU/QuEST_cpu_local.c.o"
	cd /Users/cathan/PhD/事务-课题组/Tutorial-website/2-sat/build/QuEST && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -MD -MT QuEST/CMakeFiles/QuEST.dir/src/CPU/QuEST_cpu_local.c.o -MF CMakeFiles/QuEST.dir/src/CPU/QuEST_cpu_local.c.o.d -o CMakeFiles/QuEST.dir/src/CPU/QuEST_cpu_local.c.o -c /Users/cathan/PhD/事务-课题组/Tutorial-website/2-sat/QuEST/src/CPU/QuEST_cpu_local.c

QuEST/CMakeFiles/QuEST.dir/src/CPU/QuEST_cpu_local.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/QuEST.dir/src/CPU/QuEST_cpu_local.c.i"
	cd /Users/cathan/PhD/事务-课题组/Tutorial-website/2-sat/build/QuEST && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /Users/cathan/PhD/事务-课题组/Tutorial-website/2-sat/QuEST/src/CPU/QuEST_cpu_local.c > CMakeFiles/QuEST.dir/src/CPU/QuEST_cpu_local.c.i

QuEST/CMakeFiles/QuEST.dir/src/CPU/QuEST_cpu_local.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/QuEST.dir/src/CPU/QuEST_cpu_local.c.s"
	cd /Users/cathan/PhD/事务-课题组/Tutorial-website/2-sat/build/QuEST && /Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /Users/cathan/PhD/事务-课题组/Tutorial-website/2-sat/QuEST/src/CPU/QuEST_cpu_local.c -o CMakeFiles/QuEST.dir/src/CPU/QuEST_cpu_local.c.s

# Object files for target QuEST
QuEST_OBJECTS = \
"CMakeFiles/QuEST.dir/src/QuEST.c.o" \
"CMakeFiles/QuEST.dir/src/QuEST_common.c.o" \
"CMakeFiles/QuEST.dir/src/QuEST_qasm.c.o" \
"CMakeFiles/QuEST.dir/src/QuEST_validation.c.o" \
"CMakeFiles/QuEST.dir/src/mt19937ar.c.o" \
"CMakeFiles/QuEST.dir/src/CPU/QuEST_cpu.c.o" \
"CMakeFiles/QuEST.dir/src/CPU/QuEST_cpu_local.c.o"

# External object files for target QuEST
QuEST_EXTERNAL_OBJECTS =

QuEST/libQuEST.dylib: QuEST/CMakeFiles/QuEST.dir/src/QuEST.c.o
QuEST/libQuEST.dylib: QuEST/CMakeFiles/QuEST.dir/src/QuEST_common.c.o
QuEST/libQuEST.dylib: QuEST/CMakeFiles/QuEST.dir/src/QuEST_qasm.c.o
QuEST/libQuEST.dylib: QuEST/CMakeFiles/QuEST.dir/src/QuEST_validation.c.o
QuEST/libQuEST.dylib: QuEST/CMakeFiles/QuEST.dir/src/mt19937ar.c.o
QuEST/libQuEST.dylib: QuEST/CMakeFiles/QuEST.dir/src/CPU/QuEST_cpu.c.o
QuEST/libQuEST.dylib: QuEST/CMakeFiles/QuEST.dir/src/CPU/QuEST_cpu_local.c.o
QuEST/libQuEST.dylib: QuEST/CMakeFiles/QuEST.dir/build.make
QuEST/libQuEST.dylib: QuEST/CMakeFiles/QuEST.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/cathan/PhD/事务-课题组/Tutorial-website/2-sat/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Linking C shared library libQuEST.dylib"
	cd /Users/cathan/PhD/事务-课题组/Tutorial-website/2-sat/build/QuEST && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/QuEST.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
QuEST/CMakeFiles/QuEST.dir/build: QuEST/libQuEST.dylib
.PHONY : QuEST/CMakeFiles/QuEST.dir/build

QuEST/CMakeFiles/QuEST.dir/clean:
	cd /Users/cathan/PhD/事务-课题组/Tutorial-website/2-sat/build/QuEST && $(CMAKE_COMMAND) -P CMakeFiles/QuEST.dir/cmake_clean.cmake
.PHONY : QuEST/CMakeFiles/QuEST.dir/clean

QuEST/CMakeFiles/QuEST.dir/depend:
	cd /Users/cathan/PhD/事务-课题组/Tutorial-website/2-sat/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/cathan/PhD/事务-课题组/Tutorial-website/2-sat /Users/cathan/PhD/事务-课题组/Tutorial-website/2-sat/QuEST /Users/cathan/PhD/事务-课题组/Tutorial-website/2-sat/build /Users/cathan/PhD/事务-课题组/Tutorial-website/2-sat/build/QuEST /Users/cathan/PhD/事务-课题组/Tutorial-website/2-sat/build/QuEST/CMakeFiles/QuEST.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : QuEST/CMakeFiles/QuEST.dir/depend
