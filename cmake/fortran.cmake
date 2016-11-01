

if (DEFINED USE_FORTRAN)
  set(SPLATT_NOWARN ${USE_FORTRAN})

  # Enable linking against Fortran
  enable_language(Fortran)

  # Link against a supplied Fortran library.
  if (DEFINED USER_FORTRAN_LIB)
    message(STATUS "Using user supplied Fortran library=${USER_FORTRAN_LIB}")
    set(SPLATT_LIBS ${SPLATT_LIBS} ${USER_FORTRAN_LIB})

  else()

    # Try popular ones.
    if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
      set(FORT_LIB ifcore)
    else()
      set(FORT_LIB gfortran)
    endif()
    
    # Find and use the library.
    find_library(FORTRAN_LIB ${FORT_LIB})
    set(SPLATT_LIBS ${SPLATT_LIBS} ${FORTRAN_LIB})
  endif()
endif()
