# CMake build script for the IBM ILOG Concert/CP solver for AMPL.

find_package(CPLEX)
if (CPLEX_FOUND AND CPLEX_CONCERT_FOUND AND CPLEX_CP_FOUND)
  add_definitions(-DILOGCP_NO_VERS)
  add_definitions(${CPLEX_CONCERT_DEFINITIONS} ${CPLEX_ILOCPLEX_DEFINITIONS})
  include_directories(.. ${CPLEX_ILOCPLEX_INCLUDE_DIRS}
    ${CPLEX_CP_INCLUDE_DIRS})

  add_library(amplilogcp-static concert.cc concert.h
    ilogcp.cc ilogcp.h ilogcp_date.h ../util/util-all.cc)
  target_link_libraries(amplilogcp-static
    ${CPLEX_CP_LIBRARIES} ${CPLEX_ILOCPLEX_LIBRARIES} ${CPLEX_LIBRARIES})
  add_dependencies(amplilogcp-static arith_h)

  # Link with librt for clock_gettime (Linux on i386).
  find_library(RT_LIBRARY rt)
  if (RT_LIBRARY)
    target_link_libraries(amplilogcp-static ${RT_LIBRARY})
  endif ()

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
    target_link_libraries(amplilogcp-static amplsolver-dynrt)
  else ()
    target_link_libraries(amplilogcp-static amplsolver)
  endif ()

  add_library(amplilogcp SHARED ../util/solver-c.cc)
  target_link_libraries(amplilogcp amplilogcp-static)

  add_executable(ilogcp main.cc)
  target_link_libraries(ilogcp amplilogcp-static)
  
  function(copy_cplex_dll target)
    if (CPLEX_DLL)
      get_filename_component(CPLEX_DLL_NAME "${CPLEX_DLL}" NAME)
      add_custom_command(TARGET ${target}
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${CPLEX_DLL}" "$<TARGET_FILE_DIR:${target}>/${CPLEX_DLL_NAME}")
    endif ()
  endfunction()
  copy_cplex_dll(ilogcp)
endif ()
