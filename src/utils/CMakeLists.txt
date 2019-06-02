# bar library
SET(UTILS_src ${CMAKE_SOURCE_DIR}/src/utils/dnrm2.f ${CMAKE_SOURCE_DIR}/src/utils/dpmeps.f ${CMAKE_SOURCE_DIR}/src/utils/surn01.f)
ADD_LIBRARY(utils STATIC ${UTILS_src})
set_target_properties(utils PROPERTIES LINKER_LANGUAGE Fortran)
