# CMake build script for the IBM ILOG Concert/CP solver for AMPL.

find_package(CPLEX)
if (CPLEX_FOUND AND CPLEX_CONCERT_FOUND AND CPLEX_CP_FOUND)
  add_definitions(${CPLEX_CONCERT_DEFINITIONS} ${CPLEX_ILOCPLEX_DEFINITIONS})
  include_directories(.. ${CPLEX_ILOCPLEX_INCLUDE_DIRS}
    ${CPLEX_CP_INCLUDE_DIRS})

  add_library(amplilogcp
    ilogcp.cpp ilogcp.h ilogcp_date.h ../util/util-all.cc)
  target_link_libraries(amplilogcp
    ${CPLEX_ILOCPLEX_LIBRARIES} ${CPLEX_CP_LIBRARIES} ${CPLEX_LIBRARIES})
  add_dependencies(amplilogcp arith_h)

  if (MSVC)
    # CP Optimizer only provides libraries built with DLL MSVC runtime
    # so change the compiler flags and build a compatible version of the
    # amplsolver library specifically for ilogcp.
    use_dll_runtime()
    foreach (src ${ASL_SOURCES})
      get_filename_component(path ${src} PATH)
      if (path STREQUAL "")
        set(src ../${src})
      endif ()
      set(SOURCES ${SOURCES} ${src})
    endforeach ()
    add_library(amplsolver-dynrt ${SOURCES})
    add_dependencies(amplsolver-dynrt arith_h)
    target_link_libraries(amplilogcp amplsolver-dynrt)
  else ()
    target_link_libraries(amplilogcp ${SOLVER_LIBS})
  endif ()

  add_executable(ilogcp main.cpp)
  target_link_libraries(ilogcp amplilogcp)
endif ()
