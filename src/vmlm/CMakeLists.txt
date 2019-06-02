# bar library
SET(VMLM_src ${CMAKE_SOURCE_DIR}/src/vmlm/dlmmv.f ${CMAKE_SOURCE_DIR}/src/vmlm/dvmlm.f)
ADD_LIBRARY(vmlm STATIC ${VMLM_src})
set_target_properties(vmlm PROPERTIES LINKER_LANGUAGE Fortran)
