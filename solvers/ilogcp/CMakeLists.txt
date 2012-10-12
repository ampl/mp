# CMake build script for the IBM ILOG Concert/CP driver for AMPL.

find_package(CPLEX)
if (CPLEX_FOUND AND CPLEX_CONCERT_FOUND AND CPLEX_CP_FOUND)
  add_definitions(${CPLEX_CONCERT_DEFINITIONS} ${CPLEX_ILOCPLEX_DEFINITIONS})
  include_directories(.. ${CPLEX_ILOCPLEX_INCLUDE_DIRS}
    ${CPLEX_CP_INCLUDE_DIRS})

  if (MSVC)
    # CP Optimizer only provides libraries built with DLL MSVC runtime
    # so change the compiler flags and build a compatible version of the
    # amplsolver library specifically for ilogcp.
    foreach (var
        CMAKE_C_FLAGS CMAKE_C_FLAGS_DEBUG CMAKE_C_FLAGS_RELEASE
        CMAKE_C_FLAGS_MINSIZEREL CMAKE_C_FLAGS_RELWITHDEBINFO
        CMAKE_CXX_FLAGS CMAKE_CXX_FLAGS_DEBUG CMAKE_CXX_FLAGS_RELEASE
        CMAKE_CXX_FLAGS_MINSIZEREL CMAKE_CXX_FLAGS_RELWITHDEBINFO)
      if (${var} MATCHES "${match}")
        string(REPLACE "/MT" "/MD" ${var} "${${var}}")
      endif ()
    endforeach ()
    foreach (src ${ASL_SOURCES})
      get_filename_component(path ${src} PATH)
      if (path STREQUAL "")
        set(src ../${src})
      endif ()
      set(SOURCES ${SOURCES} ${src})
    endforeach ()
    add_library(amplsolver-dynrt ${SOURCES})
    set(SOLVER_LIBS amplsolver-dynrt)
  endif ()

  add_library(amplilogcp
    build_expr.cpp
    ilogcp.cpp
    ilogcp.h
    ilogcp_date.h
    util.cpp
    util.h)
  add_dependencies(amplilogcp arith_h)

  add_executable(ilogcp main.cpp)

  target_link_libraries(ilogcp amplilogcp ${SOLVER_LIBS}
    ${CPLEX_ILOCPLEX_LIBRARIES} ${CPLEX_CP_LIBRARIES} ${CPLEX_LIBRARIES})
endif ()
