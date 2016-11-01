
include(CheckFunctionExists)

# just use MKL flag
if (INTEL_OPT)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mkl")
  add_definitions(-DSPLATT_INC_LAPACKE="mkl_lapacke.h")
  add_definitions(-DSPLATT_INC_CBLAS="mkl_cblas.h")
else()

  # Download BLAS/LAPACK
  if (DEFINED DOWNLOAD_BLAS_LAPACK)

    # Enable linking against Fortran
    enable_language(Fortran)

    message(STATUS "Downloading OpenBLAS...")
    if(${OPENMP_FOUND})
      set(BLAS_THREAD "USE_OPENMP=1")
    else()
      set(BLAS_THREAD "USE_THREADS=1")
    endif()
    execute_process(
        COMMAND ${CMAKE_SOURCE_DIR}/scripts/download-blas-lapack.sh
            ${CMAKE_BINARY_DIR}
            CC=${CMAKE_C_COMPILER}
            FC=${CMAKE_Fortran_COMPILER}
            ${BLAS_THREAD})

    set(LAPACK_LIBRARIES ${CMAKE_BINARY_DIR}/lapack/lib/libopenblas.a)
    set(SPLATT_INCLUDES ${SPLATT_INCLUDES} ${CMAKE_BINARY_DYR}/lapack/include)

    # Make sure to link against Fortran libs later.
    set(USE_FORTRAN TRUE)

    # avoid annoying warning
    set(SPLATT_NOWARN ${DOWNLOAD_BLAS_LAPACK})
  endif() # end download



  if (DEFINED USER_LAPACK_LIB)
    message("Using user supplied LAPACK=${USER_LAPACK_LIB}")
    set(LAPACK_LIBRARIES ${USER_LAPACK_LIB})
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

  # Explicitly use BLAS if user supplied. Most LAPACk libs have BLAS already...
  if (DEFINED USER_BLAS_LIB)
    message("Using user supplied BLAS=${USER_BLAS_LIB}")
    set(BLAS_LIBRARIES ${USER_BLAS_LIB})
  endif()

  #
  # Make sure we actually have C bindings
  #
  set(CMAKE_REQUIRED_LIBRARIES ${LAPACK_LIBRARIES})
  CHECK_FUNCTION_EXISTS(LAPACKE_dpotrf HAVE_LAPACKE)
  if(NOT ${HAVE_LAPACKE})
    message(FATAL_ERROR "LAPACK library does not export C bindings (LAPACKE). Run `./configure --help for assistance.")
  endif()

  set(CMAKE_REQUIRED_LIBRARIES ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
  CHECK_FUNCTION_EXISTS(cblas_dgemm HAVE_CBLAS)
  if(NOT ${HAVE_CBLAS})
    message(FATAL_ERROR "BLAS library does not export C bindings (CBLAS). Run `./configure --help for assistance.")
  endif()
  set(CMAKE_REQUIRED_LIBRARIES "")

  #
  # Account for MKL using a different header.
  #
  if("${LAPACK_LIBRARIES}" MATCHES "mkl")
    add_definitions(-DSPLATT_INC_LAPACKE="mkl_lapacke.h")
    add_definitions(-DSPLATT_INC_CBLAS="mkl_cblas.h")
  endif()

  #
  # Finally, add libraries to SPLATT_LIB.
  #
  set(SPLATT_LIBS ${SPLATT_LIBS} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
endif() # not INTEL_OPT

