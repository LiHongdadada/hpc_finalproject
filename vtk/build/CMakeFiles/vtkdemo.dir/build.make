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
CMAKE_COMMAND = /work/mae-mawb/cmake-3.20.2/cmake_install_files/bin/cmake

# The command to remove a file.
RM = /work/mae-mawb/cmake-3.20.2/cmake_install_files/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /work/mae-mawb/hpc_finalproject/vtk

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /work/mae-mawb/hpc_finalproject/vtk/build

# Include any dependencies generated for this target.
include CMakeFiles/vtkdemo.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/vtkdemo.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/vtkdemo.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/vtkdemo.dir/flags.make

CMakeFiles/vtkdemo.dir/vtk.cpp.o: CMakeFiles/vtkdemo.dir/flags.make
CMakeFiles/vtkdemo.dir/vtk.cpp.o: ../vtk.cpp
CMakeFiles/vtkdemo.dir/vtk.cpp.o: CMakeFiles/vtkdemo.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/work/mae-mawb/hpc_finalproject/vtk/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/vtkdemo.dir/vtk.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/vtkdemo.dir/vtk.cpp.o -MF CMakeFiles/vtkdemo.dir/vtk.cpp.o.d -o CMakeFiles/vtkdemo.dir/vtk.cpp.o -c /work/mae-mawb/hpc_finalproject/vtk/vtk.cpp

CMakeFiles/vtkdemo.dir/vtk.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/vtkdemo.dir/vtk.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /work/mae-mawb/hpc_finalproject/vtk/vtk.cpp > CMakeFiles/vtkdemo.dir/vtk.cpp.i

CMakeFiles/vtkdemo.dir/vtk.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/vtkdemo.dir/vtk.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /work/mae-mawb/hpc_finalproject/vtk/vtk.cpp -o CMakeFiles/vtkdemo.dir/vtk.cpp.s

# Object files for target vtkdemo
vtkdemo_OBJECTS = \
"CMakeFiles/vtkdemo.dir/vtk.cpp.o"

# External object files for target vtkdemo
vtkdemo_EXTERNAL_OBJECTS =

vtkdemo: CMakeFiles/vtkdemo.dir/vtk.cpp.o
vtkdemo: CMakeFiles/vtkdemo.dir/build.make
vtkdemo: /work/mae-mawb/VTK_build/lib/libvtkIOLegacy-8.2.so.1
vtkdemo: /work/mae-mawb/VTK_build/lib/libvtkIOXML-8.2.so.1
vtkdemo: /work/mae-mawb/VTK_build/lib/libvtkIOXMLParser-8.2.so.1
vtkdemo: /work/mae-mawb/VTK_build/lib/libvtkexpat-8.2.so.1
vtkdemo: /share/base/hdf5/1.10.4-gcc-4.8.5/lib/libhdf5.so
vtkdemo: /share/base/szip/2.1.1-gcc-4.8.5/lib/libsz.so
vtkdemo: /share/base/zlib/1.2.11-gcc-4.8.5/lib/libz.so
vtkdemo: /usr/lib64/libdl.so
vtkdemo: /usr/lib64/libm.so
vtkdemo: /work/mae-mawb/VTK_build/lib/libvtkIOCore-8.2.so.1
vtkdemo: /work/mae-mawb/VTK_build/lib/libvtkCommonExecutionModel-8.2.so.1
vtkdemo: /work/mae-mawb/VTK_build/lib/libvtkdoubleconversion-8.2.so.1
vtkdemo: /work/mae-mawb/VTK_build/lib/libvtklz4-8.2.so.1
vtkdemo: /work/mae-mawb/VTK_build/lib/libvtklzma-8.2.so.1
vtkdemo: /work/mae-mawb/VTK_build/lib/libvtkzlib-8.2.so.1
vtkdemo: /work/mae-mawb/VTK_build/lib/libvtkCommonDataModel-8.2.so.1
vtkdemo: /work/mae-mawb/VTK_build/lib/libvtkCommonSystem-8.2.so.1
vtkdemo: /work/mae-mawb/VTK_build/lib/libvtkCommonMisc-8.2.so.1
vtkdemo: /work/mae-mawb/VTK_build/lib/libvtkCommonTransforms-8.2.so.1
vtkdemo: /work/mae-mawb/VTK_build/lib/libvtkCommonMath-8.2.so.1
vtkdemo: /work/mae-mawb/VTK_build/lib/libvtkCommonCore-8.2.so.1
vtkdemo: /work/mae-mawb/VTK_build/lib/libvtksys-8.2.so.1
vtkdemo: CMakeFiles/vtkdemo.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/work/mae-mawb/hpc_finalproject/vtk/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable vtkdemo"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/vtkdemo.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/vtkdemo.dir/build: vtkdemo
.PHONY : CMakeFiles/vtkdemo.dir/build

CMakeFiles/vtkdemo.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/vtkdemo.dir/cmake_clean.cmake
.PHONY : CMakeFiles/vtkdemo.dir/clean

CMakeFiles/vtkdemo.dir/depend:
	cd /work/mae-mawb/hpc_finalproject/vtk/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /work/mae-mawb/hpc_finalproject/vtk /work/mae-mawb/hpc_finalproject/vtk /work/mae-mawb/hpc_finalproject/vtk/build /work/mae-mawb/hpc_finalproject/vtk/build /work/mae-mawb/hpc_finalproject/vtk/build/CMakeFiles/vtkdemo.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/vtkdemo.dir/depend
