# bar library
SET(BLAS_src ${CMAKE_SOURCE_DIR}/src/blas/dasum.f ${CMAKE_SOURCE_DIR}/src/blas/daxpy.f ${CMAKE_SOURCE_DIR}/src/blas/dcopy.f ${CMAKE_SOURCE_DIR}/src/blas/ddot.f ${CMAKE_SOURCE_DIR}/src/blas/dgemm.f ${CMAKE_SOURCE_DIR}/src/blas/dgemv.f 
            ${CMAKE_SOURCE_DIR}/src/blas/dger.f ${CMAKE_SOURCE_DIR}/src/blas/dnrm2.f ${CMAKE_SOURCE_DIR}/src/blas/dscal.f ${CMAKE_SOURCE_DIR}/src/blas/dswap.f ${CMAKE_SOURCE_DIR}/src/blas/dsyrk.f ${CMAKE_SOURCE_DIR}/src/blas/dtrmm.f 
            ${CMAKE_SOURCE_DIR}/src/blas/dtrmv.f ${CMAKE_SOURCE_DIR}/src/blas/dtrsm.f ${CMAKE_SOURCE_DIR}/src/blas/dtrsv.f ${CMAKE_SOURCE_DIR}/src/blas/idamax.f ${CMAKE_SOURCE_DIR}/src/blas/lsame.f ${CMAKE_SOURCE_DIR}/src/blas/xerbla.f)
ADD_LIBRARY(blas STATIC ${BLAS_src})
set_target_properties(blas PROPERTIES LINKER_LANGUAGE Fortran)
