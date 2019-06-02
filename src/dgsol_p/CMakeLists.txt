########################################
# Set up how to compile the source files
########################################
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

# Add the source files
SET(DGSOLP_src ${CMAKE_SOURCE_DIR}/src/dgsol_p/dgsol.c
)

# Define the executable in terms of the source files
ADD_EXECUTABLE(dgsol_p ${DGSOLP_src})
set_target_properties(dgsol_s PROPERTIES LINKER_LANGUAGE C)

#####################################################
# Add the needed libraries and special compiler flags
#####################################################

# This links foo to the bar library
TARGET_LINK_LIBRARIES(dgsol_p  dgvmlm vmlm csrch  blas utils ${MPI_LIBRARIES})
link_libraries(dgsol_p gfortran)

# Uncomment if you need to link to BLAS and LAPACK
#TARGET_LINK_LIBRARIES(${FOOEXE} ${BLAS_LIBRARIES}
#                                ${LAPACK_LIBRARIES}
#                                ${CMAKE_THREAD_LIBS_INIT})

# Uncomment if you have parallization
#IF(USE_OPENMP)
#    SET_TARGET_PROPERTIES(${FOOEXE} PROPERTIES
#                          COMPILE_FLAGS "${OpenMP_Fortran_FLAGS}"
#                          LINK_FLAGS "${OpenMP_Fortran_FLAGS}")
#ELSEIF(USE_MPI)
#    SET_TARGET_PROPERTIES(${FOOEXE} PROPERTIES
#                          COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS}"
#                          LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}")
#    INCLUDE_DIRECTORIES(${MPI_Fortran_INCLUDE_PATH})
#    TARGET_LINK_LIBRARIES(${FOOEXE} ${MPI_Fortran_LIBRARIES})
#ENDIF(USE_OPENMP)

#####################################
# Tell how to install this executable
#####################################

