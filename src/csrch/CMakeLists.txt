# bar library
SET(CSRCH_src ${CMAKE_SOURCE_DIR}/src/csrch/dcsrch.f ${CMAKE_SOURCE_DIR}/src/csrch/dcstep.f)
ADD_LIBRARY(csrch STATIC ${CSRCH_src})
set_target_properties(csrch PROPERTIES LINKER_LANGUAGE Fortran)
