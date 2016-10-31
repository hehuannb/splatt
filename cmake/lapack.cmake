
include(CheckFunctionExists)

# just use MKL flag
if (INTEL_OPT)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mkl")
  add_definitions(-DSPLATT_INC_LAPACKE="mkl_lapacke.h")
else()

  # Download BLAS/LAPACK
  if (DEFINED DOWNLOAD_BLAS_LAPACK)

    # Enable linking against Fortran
    enable_language(Fortran)

    message(WARNING "Downloading OpenBLAS.")
    if(${OPENMP_FOUND})
      execute_process(COMMAND ${CMAKE_SOURCE_DIR}/scripts/download-blas-lapack.sh ${CMAKE_BINARY_DIR} CC=${CMAKE_C_COMPILER} FC=${CMAKE_Fortran_COMPILER} USE_OPENMP=1)
    else()
      execute_process(COMMAND ${CMAKE_SOURCE_DIR}/scripts/download-blas-lapack.sh ${CMAKE_BINARY_DIR} CC=${CMAKE_C_COMPILER} FC=${CMAKE_Fortran_COMPILER} USE_THREADS=1)
    endif()

    set(USER_LAPACK_LIB ${CMAKE_BINARY_DIR}/lapack/lib/libopenblas.a)
    set(USER_BLAS_LIB ${CMAKE_BINARY_DIR}/lapack/lib/libopenblas.a)
    set(SPLATT_INCLUDES ${SPLATT_INCLUDES} ${CMAKE_BINARY_DYR}/lapack/include)

    # Link against generic BLAS and a Fortran library.
    # TODO: Is there a better way to do this? The Fortran library must be added
    # AFTER BLAS/LAPACK.
    if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
      set(USER_BLAS_LIB ${USER_BLAS_LIB} ifcore)
    else()
      set(USER_BLAS_LIB ${USER_BLAS_LIB} gfortran)
    endif()

    # avoid annoying warning
    set(SPLATT_NOWARN ${DOWNLOAD_BLAS_LAPACK})
  endif()

  if (DEFINED USER_LAPACK_LIB)
    message("Using user supplied LAPACK=${USER_LAPACK_LIB}")
    set(SPLATT_LIBS ${SPLATT_LIBS} ${USER_LAPACK_LIB})
  # auto find LAPACK
  else()
    find_package(LAPACK)
    set(SPLATT_LIBS ${SPLATT_LIBS} ${LAPACK_LIBRARIES})
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}")
    if(${LAPACK_FOUND})
      message("FOUND LAPACK LIBS: " ${LAPACK_LIBRARIES})
    else()
      message(FATAL_ERROR "Could not find LAPACK library. Run `./configure --help`  for assistance.")
    endif()
  endif()

  if (DEFINED USER_BLAS_LIB)
    message("Using user supplied BLAS=${USER_BLAS_LIB}")
    set(SPLATT_LIBS ${SPLATT_LIBS} ${USER_BLAS_LIB})
  # auto find BLAS
  else()
    find_package(BLAS)
    set(SPLATT_LIBS ${SPLATT_LIBS} ${BLAS_LIBRARIES})
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${BLAS_LINKER_FLAGS}")
    if(${BLAS_FOUND})
      message("FOUND BLAS LIBS: " ${BLAS_LIBRARIES})
    else()
      message(FATAL_ERROR "Could not find BLAS library. Run `./configure --help`  for assistance.")
    endif()
  endif()



  # Make sure we actually have C bindings
  set(CMAKE_REQUIRED_LIBRARIES ${SPLATT_LIBS})
  CHECK_FUNCTION_EXISTS(LAPACKE_dpotrf HAVE_LAPACKE)
  if(NOT ${HAVE_LAPACKE})
    message(FATAL_ERROR "LAPACK library does not export C bindings (LAPACKE). Run `./configure --help for assistance.")
  endif()

  CHECK_FUNCTION_EXISTS(cblas_dgemm HAVE_CBLAS)
  if(NOT ${HAVE_CBLAS})
    message(FATAL_ERROR "BLAS library does not export C bindings (CBLAS). Run `./configure --help for assistance.")
  endif()
  set(CMAKE_REQUIRED_LIBRARIES "")

endif() # not INTEL_OPT

