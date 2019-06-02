# CMake project file for FOO

##################################################
# Define the project and the depencies that it has
##################################################

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.5)
#set(CMAKE_VERBOSE_MAKEFILE ON)
PROJECT(DGSOL Fortran CXX C)

# Set the FOO version
SET(VERSION 1.2.3)

# Add our local modlues to the module path
SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# Uncomment if it is required that Fortran 90 is supported
#IF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)
#    MESSAGE(FATAL_ERROR "Fortran compiler does not support F90")
#ENDIF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)

# Set some options the user may choose
# Uncomment the below if you want the user to choose a parallelization library
#OPTION(USE_MPI "Use the MPI library for parallelization" OFF)
#OPTION(USE_OPENMP "Use OpenMP for parallelization" OFF)

# This INCLUDE statement executes code that sets the compile flags for DEBUG,
# RELEASE, and TESTING.  You should  review this file and make sure the flags 
# are to your liking.
INCLUDE(${CMAKE_MODULE_PATH}/SetFortranFlags.cmake) 
# Locate and set parallelization libraries.  There are some CMake peculiarities
# taken care of here, such as the fact that the FindOpenMP routine doesn't know
# about Fortran.
#INCLUDE(${CMAKE_MODULE_PATH}/SetParallelizationLibrary.cmake)
# Setup the LAPACK libraries.  This also takes care of peculiarities, such as
# the fact the searching for MKL requires a C compiler, and that the results
# are not stored in the cache. 
#INCLUDE(${CMAKE_MODULE_PATH}/SetUpLAPACK.cmake)

# There is an error in CMAKE with this flag for pgf90.  Unset it
GET_FILENAME_COMPONENT(FCNAME ${CMAKE_Fortran_COMPILER} NAME)
IF(FCNAME STREQUAL "pgf90")
    UNSET(CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS)
ENDIF(FCNAME STREQUAL "pgf90")

############################################################
# Define the actual files and folders that make up the build
############################################################

# Define the executable name
#SET(FOOEXE foo)

# Define some directories
SET(SRC ${CMAKE_SOURCE_DIR}/src)

# Have the .mod files placed in the lib folder
SET(CMAKE_Fortran_MODULE_DIRECTORY ${LIB})

# The source for the BAR library and have it placed in the lib folder
ADD_SUBDIRECTORY(${CMAKE_SOURCE_DIR}/src/blas)
ADD_SUBDIRECTORY(${CMAKE_SOURCE_DIR}/src/csrch)
ADD_SUBDIRECTORY(${CMAKE_SOURCE_DIR}/src/dgvmlm)
ADD_SUBDIRECTORY(${CMAKE_SOURCE_DIR}/src/utils)
ADD_SUBDIRECTORY(${CMAKE_SOURCE_DIR}/src/vmlm)
ADD_SUBDIRECTORY(${CMAKE_SOURCE_DIR}/src/dgsol_s)
ADD_SUBDIRECTORY(${CMAKE_SOURCE_DIR}/src/dgsol_p)

ADD_EXECUTABLE(dgsol src/dgsol_j/dgsol_seq_lib.c src/dgsol_j/dgsolj_seq_test.cpp)
set_target_properties(dgsol PROPERTIES LINKER_LANGUAGE C)
TARGET_LINK_LIBRARIES(dgsol  dgvmlm vmlm csrch   blas utils)
link_libraries(dgsol gfortran m)

# The source for the FOO binary and have it placed in the bin folder
#ADD_SUBDIRECTORY(${SRCFOO} ${BIN})

