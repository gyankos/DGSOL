# bar library
SET(DGVMLM_src ${CMAKE_SOURCE_DIR}/src/dgvmlm/dgopt.f ${CMAKE_SOURCE_DIR}/src/dgvmlm/dginitial.f ${CMAKE_SOURCE_DIR}/src/dgvmlm/dgfun.f ${CMAKE_SOURCE_DIR}/src/dgvmlm/dgerr.f ${CMAKE_SOURCE_DIR}/src/dgvmlm/dgsel.f)
ADD_LIBRARY(dgvmlm STATIC ${DGVMLM_src})
set_target_properties(dgvmlm PROPERTIES LINKER_LANGUAGE Fortran)
