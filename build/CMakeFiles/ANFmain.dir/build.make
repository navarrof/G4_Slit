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
CMAKE_SOURCE_DIR = /mnt/c/Users/navarrof/Araceli_NF/2021/October/SlitThermalSim/SLIT_G4

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /mnt/c/Users/navarrof/Araceli_NF/2021/October/SlitThermalSim/SLIT_G4/build

# Include any dependencies generated for this target.
include CMakeFiles/ANFmain.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/ANFmain.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/ANFmain.dir/flags.make

CMakeFiles/ANFmain.dir/ANF_Main.cc.o: CMakeFiles/ANFmain.dir/flags.make
CMakeFiles/ANFmain.dir/ANF_Main.cc.o: ../ANF_Main.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/c/Users/navarrof/Araceli_NF/2021/October/SlitThermalSim/SLIT_G4/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/ANFmain.dir/ANF_Main.cc.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/ANFmain.dir/ANF_Main.cc.o -c /mnt/c/Users/navarrof/Araceli_NF/2021/October/SlitThermalSim/SLIT_G4/ANF_Main.cc

CMakeFiles/ANFmain.dir/ANF_Main.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ANFmain.dir/ANF_Main.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/c/Users/navarrof/Araceli_NF/2021/October/SlitThermalSim/SLIT_G4/ANF_Main.cc > CMakeFiles/ANFmain.dir/ANF_Main.cc.i

CMakeFiles/ANFmain.dir/ANF_Main.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ANFmain.dir/ANF_Main.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/c/Users/navarrof/Araceli_NF/2021/October/SlitThermalSim/SLIT_G4/ANF_Main.cc -o CMakeFiles/ANFmain.dir/ANF_Main.cc.s

CMakeFiles/ANFmain.dir/ANF_Main.cc.o.requires:

.PHONY : CMakeFiles/ANFmain.dir/ANF_Main.cc.o.requires

CMakeFiles/ANFmain.dir/ANF_Main.cc.o.provides: CMakeFiles/ANFmain.dir/ANF_Main.cc.o.requires
	$(MAKE) -f CMakeFiles/ANFmain.dir/build.make CMakeFiles/ANFmain.dir/ANF_Main.cc.o.provides.build
.PHONY : CMakeFiles/ANFmain.dir/ANF_Main.cc.o.provides

CMakeFiles/ANFmain.dir/ANF_Main.cc.o.provides.build: CMakeFiles/ANFmain.dir/ANF_Main.cc.o


CMakeFiles/ANFmain.dir/src/ANF_ActionInitialization.cc.o: CMakeFiles/ANFmain.dir/flags.make
CMakeFiles/ANFmain.dir/src/ANF_ActionInitialization.cc.o: ../src/ANF_ActionInitialization.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/c/Users/navarrof/Araceli_NF/2021/October/SlitThermalSim/SLIT_G4/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/ANFmain.dir/src/ANF_ActionInitialization.cc.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/ANFmain.dir/src/ANF_ActionInitialization.cc.o -c /mnt/c/Users/navarrof/Araceli_NF/2021/October/SlitThermalSim/SLIT_G4/src/ANF_ActionInitialization.cc

CMakeFiles/ANFmain.dir/src/ANF_ActionInitialization.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ANFmain.dir/src/ANF_ActionInitialization.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/c/Users/navarrof/Araceli_NF/2021/October/SlitThermalSim/SLIT_G4/src/ANF_ActionInitialization.cc > CMakeFiles/ANFmain.dir/src/ANF_ActionInitialization.cc.i

CMakeFiles/ANFmain.dir/src/ANF_ActionInitialization.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ANFmain.dir/src/ANF_ActionInitialization.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/c/Users/navarrof/Araceli_NF/2021/October/SlitThermalSim/SLIT_G4/src/ANF_ActionInitialization.cc -o CMakeFiles/ANFmain.dir/src/ANF_ActionInitialization.cc.s

CMakeFiles/ANFmain.dir/src/ANF_ActionInitialization.cc.o.requires:

.PHONY : CMakeFiles/ANFmain.dir/src/ANF_ActionInitialization.cc.o.requires

CMakeFiles/ANFmain.dir/src/ANF_ActionInitialization.cc.o.provides: CMakeFiles/ANFmain.dir/src/ANF_ActionInitialization.cc.o.requires
	$(MAKE) -f CMakeFiles/ANFmain.dir/build.make CMakeFiles/ANFmain.dir/src/ANF_ActionInitialization.cc.o.provides.build
.PHONY : CMakeFiles/ANFmain.dir/src/ANF_ActionInitialization.cc.o.provides

CMakeFiles/ANFmain.dir/src/ANF_ActionInitialization.cc.o.provides.build: CMakeFiles/ANFmain.dir/src/ANF_ActionInitialization.cc.o


CMakeFiles/ANFmain.dir/src/ANF_DetectorConstruction.cc.o: CMakeFiles/ANFmain.dir/flags.make
CMakeFiles/ANFmain.dir/src/ANF_DetectorConstruction.cc.o: ../src/ANF_DetectorConstruction.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/c/Users/navarrof/Araceli_NF/2021/October/SlitThermalSim/SLIT_G4/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/ANFmain.dir/src/ANF_DetectorConstruction.cc.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/ANFmain.dir/src/ANF_DetectorConstruction.cc.o -c /mnt/c/Users/navarrof/Araceli_NF/2021/October/SlitThermalSim/SLIT_G4/src/ANF_DetectorConstruction.cc

CMakeFiles/ANFmain.dir/src/ANF_DetectorConstruction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ANFmain.dir/src/ANF_DetectorConstruction.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/c/Users/navarrof/Araceli_NF/2021/October/SlitThermalSim/SLIT_G4/src/ANF_DetectorConstruction.cc > CMakeFiles/ANFmain.dir/src/ANF_DetectorConstruction.cc.i

CMakeFiles/ANFmain.dir/src/ANF_DetectorConstruction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ANFmain.dir/src/ANF_DetectorConstruction.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/c/Users/navarrof/Araceli_NF/2021/October/SlitThermalSim/SLIT_G4/src/ANF_DetectorConstruction.cc -o CMakeFiles/ANFmain.dir/src/ANF_DetectorConstruction.cc.s

CMakeFiles/ANFmain.dir/src/ANF_DetectorConstruction.cc.o.requires:

.PHONY : CMakeFiles/ANFmain.dir/src/ANF_DetectorConstruction.cc.o.requires

CMakeFiles/ANFmain.dir/src/ANF_DetectorConstruction.cc.o.provides: CMakeFiles/ANFmain.dir/src/ANF_DetectorConstruction.cc.o.requires
	$(MAKE) -f CMakeFiles/ANFmain.dir/build.make CMakeFiles/ANFmain.dir/src/ANF_DetectorConstruction.cc.o.provides.build
.PHONY : CMakeFiles/ANFmain.dir/src/ANF_DetectorConstruction.cc.o.provides

CMakeFiles/ANFmain.dir/src/ANF_DetectorConstruction.cc.o.provides.build: CMakeFiles/ANFmain.dir/src/ANF_DetectorConstruction.cc.o


CMakeFiles/ANFmain.dir/src/ANF_EventAction.cc.o: CMakeFiles/ANFmain.dir/flags.make
CMakeFiles/ANFmain.dir/src/ANF_EventAction.cc.o: ../src/ANF_EventAction.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/c/Users/navarrof/Araceli_NF/2021/October/SlitThermalSim/SLIT_G4/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/ANFmain.dir/src/ANF_EventAction.cc.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/ANFmain.dir/src/ANF_EventAction.cc.o -c /mnt/c/Users/navarrof/Araceli_NF/2021/October/SlitThermalSim/SLIT_G4/src/ANF_EventAction.cc

CMakeFiles/ANFmain.dir/src/ANF_EventAction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ANFmain.dir/src/ANF_EventAction.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/c/Users/navarrof/Araceli_NF/2021/October/SlitThermalSim/SLIT_G4/src/ANF_EventAction.cc > CMakeFiles/ANFmain.dir/src/ANF_EventAction.cc.i

CMakeFiles/ANFmain.dir/src/ANF_EventAction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ANFmain.dir/src/ANF_EventAction.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/c/Users/navarrof/Araceli_NF/2021/October/SlitThermalSim/SLIT_G4/src/ANF_EventAction.cc -o CMakeFiles/ANFmain.dir/src/ANF_EventAction.cc.s

CMakeFiles/ANFmain.dir/src/ANF_EventAction.cc.o.requires:

.PHONY : CMakeFiles/ANFmain.dir/src/ANF_EventAction.cc.o.requires

CMakeFiles/ANFmain.dir/src/ANF_EventAction.cc.o.provides: CMakeFiles/ANFmain.dir/src/ANF_EventAction.cc.o.requires
	$(MAKE) -f CMakeFiles/ANFmain.dir/build.make CMakeFiles/ANFmain.dir/src/ANF_EventAction.cc.o.provides.build
.PHONY : CMakeFiles/ANFmain.dir/src/ANF_EventAction.cc.o.provides

CMakeFiles/ANFmain.dir/src/ANF_EventAction.cc.o.provides.build: CMakeFiles/ANFmain.dir/src/ANF_EventAction.cc.o


CMakeFiles/ANFmain.dir/src/ANF_PrimaryGeneratorAction.cc.o: CMakeFiles/ANFmain.dir/flags.make
CMakeFiles/ANFmain.dir/src/ANF_PrimaryGeneratorAction.cc.o: ../src/ANF_PrimaryGeneratorAction.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/c/Users/navarrof/Araceli_NF/2021/October/SlitThermalSim/SLIT_G4/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/ANFmain.dir/src/ANF_PrimaryGeneratorAction.cc.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/ANFmain.dir/src/ANF_PrimaryGeneratorAction.cc.o -c /mnt/c/Users/navarrof/Araceli_NF/2021/October/SlitThermalSim/SLIT_G4/src/ANF_PrimaryGeneratorAction.cc

CMakeFiles/ANFmain.dir/src/ANF_PrimaryGeneratorAction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ANFmain.dir/src/ANF_PrimaryGeneratorAction.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/c/Users/navarrof/Araceli_NF/2021/October/SlitThermalSim/SLIT_G4/src/ANF_PrimaryGeneratorAction.cc > CMakeFiles/ANFmain.dir/src/ANF_PrimaryGeneratorAction.cc.i

CMakeFiles/ANFmain.dir/src/ANF_PrimaryGeneratorAction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ANFmain.dir/src/ANF_PrimaryGeneratorAction.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/c/Users/navarrof/Araceli_NF/2021/October/SlitThermalSim/SLIT_G4/src/ANF_PrimaryGeneratorAction.cc -o CMakeFiles/ANFmain.dir/src/ANF_PrimaryGeneratorAction.cc.s

CMakeFiles/ANFmain.dir/src/ANF_PrimaryGeneratorAction.cc.o.requires:

.PHONY : CMakeFiles/ANFmain.dir/src/ANF_PrimaryGeneratorAction.cc.o.requires

CMakeFiles/ANFmain.dir/src/ANF_PrimaryGeneratorAction.cc.o.provides: CMakeFiles/ANFmain.dir/src/ANF_PrimaryGeneratorAction.cc.o.requires
	$(MAKE) -f CMakeFiles/ANFmain.dir/build.make CMakeFiles/ANFmain.dir/src/ANF_PrimaryGeneratorAction.cc.o.provides.build
.PHONY : CMakeFiles/ANFmain.dir/src/ANF_PrimaryGeneratorAction.cc.o.provides

CMakeFiles/ANFmain.dir/src/ANF_PrimaryGeneratorAction.cc.o.provides.build: CMakeFiles/ANFmain.dir/src/ANF_PrimaryGeneratorAction.cc.o


CMakeFiles/ANFmain.dir/src/ANF_RunAction.cc.o: CMakeFiles/ANFmain.dir/flags.make
CMakeFiles/ANFmain.dir/src/ANF_RunAction.cc.o: ../src/ANF_RunAction.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/c/Users/navarrof/Araceli_NF/2021/October/SlitThermalSim/SLIT_G4/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/ANFmain.dir/src/ANF_RunAction.cc.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/ANFmain.dir/src/ANF_RunAction.cc.o -c /mnt/c/Users/navarrof/Araceli_NF/2021/October/SlitThermalSim/SLIT_G4/src/ANF_RunAction.cc

CMakeFiles/ANFmain.dir/src/ANF_RunAction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ANFmain.dir/src/ANF_RunAction.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/c/Users/navarrof/Araceli_NF/2021/October/SlitThermalSim/SLIT_G4/src/ANF_RunAction.cc > CMakeFiles/ANFmain.dir/src/ANF_RunAction.cc.i

CMakeFiles/ANFmain.dir/src/ANF_RunAction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ANFmain.dir/src/ANF_RunAction.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/c/Users/navarrof/Araceli_NF/2021/October/SlitThermalSim/SLIT_G4/src/ANF_RunAction.cc -o CMakeFiles/ANFmain.dir/src/ANF_RunAction.cc.s

CMakeFiles/ANFmain.dir/src/ANF_RunAction.cc.o.requires:

.PHONY : CMakeFiles/ANFmain.dir/src/ANF_RunAction.cc.o.requires

CMakeFiles/ANFmain.dir/src/ANF_RunAction.cc.o.provides: CMakeFiles/ANFmain.dir/src/ANF_RunAction.cc.o.requires
	$(MAKE) -f CMakeFiles/ANFmain.dir/build.make CMakeFiles/ANFmain.dir/src/ANF_RunAction.cc.o.provides.build
.PHONY : CMakeFiles/ANFmain.dir/src/ANF_RunAction.cc.o.provides

CMakeFiles/ANFmain.dir/src/ANF_RunAction.cc.o.provides.build: CMakeFiles/ANFmain.dir/src/ANF_RunAction.cc.o


CMakeFiles/ANFmain.dir/src/ANF_SteppingAction.cc.o: CMakeFiles/ANFmain.dir/flags.make
CMakeFiles/ANFmain.dir/src/ANF_SteppingAction.cc.o: ../src/ANF_SteppingAction.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/c/Users/navarrof/Araceli_NF/2021/October/SlitThermalSim/SLIT_G4/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object CMakeFiles/ANFmain.dir/src/ANF_SteppingAction.cc.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/ANFmain.dir/src/ANF_SteppingAction.cc.o -c /mnt/c/Users/navarrof/Araceli_NF/2021/October/SlitThermalSim/SLIT_G4/src/ANF_SteppingAction.cc

CMakeFiles/ANFmain.dir/src/ANF_SteppingAction.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ANFmain.dir/src/ANF_SteppingAction.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/c/Users/navarrof/Araceli_NF/2021/October/SlitThermalSim/SLIT_G4/src/ANF_SteppingAction.cc > CMakeFiles/ANFmain.dir/src/ANF_SteppingAction.cc.i

CMakeFiles/ANFmain.dir/src/ANF_SteppingAction.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ANFmain.dir/src/ANF_SteppingAction.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/c/Users/navarrof/Araceli_NF/2021/October/SlitThermalSim/SLIT_G4/src/ANF_SteppingAction.cc -o CMakeFiles/ANFmain.dir/src/ANF_SteppingAction.cc.s

CMakeFiles/ANFmain.dir/src/ANF_SteppingAction.cc.o.requires:

.PHONY : CMakeFiles/ANFmain.dir/src/ANF_SteppingAction.cc.o.requires

CMakeFiles/ANFmain.dir/src/ANF_SteppingAction.cc.o.provides: CMakeFiles/ANFmain.dir/src/ANF_SteppingAction.cc.o.requires
	$(MAKE) -f CMakeFiles/ANFmain.dir/build.make CMakeFiles/ANFmain.dir/src/ANF_SteppingAction.cc.o.provides.build
.PHONY : CMakeFiles/ANFmain.dir/src/ANF_SteppingAction.cc.o.provides

CMakeFiles/ANFmain.dir/src/ANF_SteppingAction.cc.o.provides.build: CMakeFiles/ANFmain.dir/src/ANF_SteppingAction.cc.o


# Object files for target ANFmain
ANFmain_OBJECTS = \
"CMakeFiles/ANFmain.dir/ANF_Main.cc.o" \
"CMakeFiles/ANFmain.dir/src/ANF_ActionInitialization.cc.o" \
"CMakeFiles/ANFmain.dir/src/ANF_DetectorConstruction.cc.o" \
"CMakeFiles/ANFmain.dir/src/ANF_EventAction.cc.o" \
"CMakeFiles/ANFmain.dir/src/ANF_PrimaryGeneratorAction.cc.o" \
"CMakeFiles/ANFmain.dir/src/ANF_RunAction.cc.o" \
"CMakeFiles/ANFmain.dir/src/ANF_SteppingAction.cc.o"

# External object files for target ANFmain
ANFmain_EXTERNAL_OBJECTS =

ANFmain: CMakeFiles/ANFmain.dir/ANF_Main.cc.o
ANFmain: CMakeFiles/ANFmain.dir/src/ANF_ActionInitialization.cc.o
ANFmain: CMakeFiles/ANFmain.dir/src/ANF_DetectorConstruction.cc.o
ANFmain: CMakeFiles/ANFmain.dir/src/ANF_EventAction.cc.o
ANFmain: CMakeFiles/ANFmain.dir/src/ANF_PrimaryGeneratorAction.cc.o
ANFmain: CMakeFiles/ANFmain.dir/src/ANF_RunAction.cc.o
ANFmain: CMakeFiles/ANFmain.dir/src/ANF_SteppingAction.cc.o
ANFmain: CMakeFiles/ANFmain.dir/build.make
ANFmain: /mnt/c/Users/navarrof/Araceli_NF/Programs/Geant4ANF_10_7/geant4.10.07.p02-Install/lib/libG4Tree.so
ANFmain: /mnt/c/Users/navarrof/Araceli_NF/Programs/Geant4ANF_10_7/geant4.10.07.p02-Install/lib/libG4FR.so
ANFmain: /mnt/c/Users/navarrof/Araceli_NF/Programs/Geant4ANF_10_7/geant4.10.07.p02-Install/lib/libG4GMocren.so
ANFmain: /mnt/c/Users/navarrof/Araceli_NF/Programs/Geant4ANF_10_7/geant4.10.07.p02-Install/lib/libG4visHepRep.so
ANFmain: /mnt/c/Users/navarrof/Araceli_NF/Programs/Geant4ANF_10_7/geant4.10.07.p02-Install/lib/libG4RayTracer.so
ANFmain: /mnt/c/Users/navarrof/Araceli_NF/Programs/Geant4ANF_10_7/geant4.10.07.p02-Install/lib/libG4VRML.so
ANFmain: /mnt/c/Users/navarrof/Araceli_NF/Programs/Geant4ANF_10_7/geant4.10.07.p02-Install/lib/libG4OpenGL.so
ANFmain: /mnt/c/Users/navarrof/Araceli_NF/Programs/Geant4ANF_10_7/geant4.10.07.p02-Install/lib/libG4gl2ps.so
ANFmain: /mnt/c/Users/navarrof/Araceli_NF/Programs/Geant4ANF_10_7/geant4.10.07.p02-Install/lib/libG4vis_management.so
ANFmain: /mnt/c/Users/navarrof/Araceli_NF/Programs/Geant4ANF_10_7/geant4.10.07.p02-Install/lib/libG4modeling.so
ANFmain: /mnt/c/Users/navarrof/Araceli_NF/Programs/Geant4ANF_10_7/geant4.10.07.p02-Install/lib/libG4interfaces.so
ANFmain: /mnt/c/Users/navarrof/Araceli_NF/Programs/Geant4ANF_10_7/geant4.10.07.p02-Install/lib/libG4persistency.so
ANFmain: /mnt/c/Users/navarrof/Araceli_NF/Programs/Geant4ANF_10_7/geant4.10.07.p02-Install/lib/libG4error_propagation.so
ANFmain: /mnt/c/Users/navarrof/Araceli_NF/Programs/Geant4ANF_10_7/geant4.10.07.p02-Install/lib/libG4readout.so
ANFmain: /mnt/c/Users/navarrof/Araceli_NF/Programs/Geant4ANF_10_7/geant4.10.07.p02-Install/lib/libG4physicslists.so
ANFmain: /mnt/c/Users/navarrof/Araceli_NF/Programs/Geant4ANF_10_7/geant4.10.07.p02-Install/lib/libG4tasking.so
ANFmain: /mnt/c/Users/navarrof/Araceli_NF/Programs/Geant4ANF_10_7/geant4.10.07.p02-Install/lib/libG4parmodels.so
ANFmain: /usr/lib/x86_64-linux-gnu/libGL.so
ANFmain: /usr/lib/x86_64-linux-gnu/libQt5OpenGL.so.5.9.5
ANFmain: /usr/lib/x86_64-linux-gnu/libQt5PrintSupport.so.5.9.5
ANFmain: /usr/lib/x86_64-linux-gnu/libQt5Widgets.so.5.9.5
ANFmain: /usr/lib/x86_64-linux-gnu/libQt5Gui.so.5.9.5
ANFmain: /usr/lib/x86_64-linux-gnu/libQt5Core.so.5.9.5
ANFmain: /usr/lib/x86_64-linux-gnu/libXmu.so
ANFmain: /usr/lib/x86_64-linux-gnu/libXext.so
ANFmain: /usr/lib/x86_64-linux-gnu/libXt.so
ANFmain: /usr/lib/x86_64-linux-gnu/libICE.so
ANFmain: /usr/lib/x86_64-linux-gnu/libSM.so
ANFmain: /usr/lib/x86_64-linux-gnu/libX11.so
ANFmain: /usr/lib/x86_64-linux-gnu/libxerces-c.so
ANFmain: /mnt/c/Users/navarrof/Araceli_NF/Programs/Geant4ANF_10_7/geant4.10.07.p02-Install/lib/libG4run.so
ANFmain: /mnt/c/Users/navarrof/Araceli_NF/Programs/Geant4ANF_10_7/geant4.10.07.p02-Install/lib/libG4event.so
ANFmain: /mnt/c/Users/navarrof/Araceli_NF/Programs/Geant4ANF_10_7/geant4.10.07.p02-Install/lib/libG4tracking.so
ANFmain: /mnt/c/Users/navarrof/Araceli_NF/Programs/Geant4ANF_10_7/geant4.10.07.p02-Install/lib/libG4processes.so
ANFmain: /mnt/c/Users/navarrof/Araceli_NF/Programs/Geant4ANF_10_7/geant4.10.07.p02-Install/lib/libG4analysis.so
ANFmain: /usr/lib/x86_64-linux-gnu/libexpat.so
ANFmain: /mnt/c/Users/navarrof/Araceli_NF/Programs/Geant4ANF_10_7/geant4.10.07.p02-Install/lib/libG4digits_hits.so
ANFmain: /mnt/c/Users/navarrof/Araceli_NF/Programs/Geant4ANF_10_7/geant4.10.07.p02-Install/lib/libG4track.so
ANFmain: /mnt/c/Users/navarrof/Araceli_NF/Programs/Geant4ANF_10_7/geant4.10.07.p02-Install/lib/libG4particles.so
ANFmain: /mnt/c/Users/navarrof/Araceli_NF/Programs/Geant4ANF_10_7/geant4.10.07.p02-Install/lib/libG4geometry.so
ANFmain: /mnt/c/Users/navarrof/Araceli_NF/Programs/Geant4ANF_10_7/geant4.10.07.p02-Install/lib/libG4materials.so
ANFmain: /mnt/c/Users/navarrof/Araceli_NF/Programs/Geant4ANF_10_7/geant4.10.07.p02-Install/lib/libG4zlib.so
ANFmain: /mnt/c/Users/navarrof/Araceli_NF/Programs/Geant4ANF_10_7/geant4.10.07.p02-Install/lib/libG4graphics_reps.so
ANFmain: /mnt/c/Users/navarrof/Araceli_NF/Programs/Geant4ANF_10_7/geant4.10.07.p02-Install/lib/libG4intercoms.so
ANFmain: /mnt/c/Users/navarrof/Araceli_NF/Programs/Geant4ANF_10_7/geant4.10.07.p02-Install/lib/libG4global.so
ANFmain: /mnt/c/Users/navarrof/Araceli_NF/Programs/Geant4ANF_10_7/geant4.10.07.p02-Install/lib/libG4clhep.so
ANFmain: /mnt/c/Users/navarrof/Araceli_NF/Programs/Geant4ANF_10_7/geant4.10.07.p02-Install/lib/libG4ptl.so.0.0.2
ANFmain: CMakeFiles/ANFmain.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/mnt/c/Users/navarrof/Araceli_NF/2021/October/SlitThermalSim/SLIT_G4/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Linking CXX executable ANFmain"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/ANFmain.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/ANFmain.dir/build: ANFmain

.PHONY : CMakeFiles/ANFmain.dir/build

CMakeFiles/ANFmain.dir/requires: CMakeFiles/ANFmain.dir/ANF_Main.cc.o.requires
CMakeFiles/ANFmain.dir/requires: CMakeFiles/ANFmain.dir/src/ANF_ActionInitialization.cc.o.requires
CMakeFiles/ANFmain.dir/requires: CMakeFiles/ANFmain.dir/src/ANF_DetectorConstruction.cc.o.requires
CMakeFiles/ANFmain.dir/requires: CMakeFiles/ANFmain.dir/src/ANF_EventAction.cc.o.requires
CMakeFiles/ANFmain.dir/requires: CMakeFiles/ANFmain.dir/src/ANF_PrimaryGeneratorAction.cc.o.requires
CMakeFiles/ANFmain.dir/requires: CMakeFiles/ANFmain.dir/src/ANF_RunAction.cc.o.requires
CMakeFiles/ANFmain.dir/requires: CMakeFiles/ANFmain.dir/src/ANF_SteppingAction.cc.o.requires

.PHONY : CMakeFiles/ANFmain.dir/requires

CMakeFiles/ANFmain.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/ANFmain.dir/cmake_clean.cmake
.PHONY : CMakeFiles/ANFmain.dir/clean

CMakeFiles/ANFmain.dir/depend:
	cd /mnt/c/Users/navarrof/Araceli_NF/2021/October/SlitThermalSim/SLIT_G4/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /mnt/c/Users/navarrof/Araceli_NF/2021/October/SlitThermalSim/SLIT_G4 /mnt/c/Users/navarrof/Araceli_NF/2021/October/SlitThermalSim/SLIT_G4 /mnt/c/Users/navarrof/Araceli_NF/2021/October/SlitThermalSim/SLIT_G4/build /mnt/c/Users/navarrof/Araceli_NF/2021/October/SlitThermalSim/SLIT_G4/build /mnt/c/Users/navarrof/Araceli_NF/2021/October/SlitThermalSim/SLIT_G4/build/CMakeFiles/ANFmain.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/ANFmain.dir/depend
