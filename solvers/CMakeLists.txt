# Assigns the value of input_var to output_var if input_var is defined.
# Otherwise assigns the values of ARGN to output_var.
function(get_var output_var input_var)
  if (${input_var})
    set(${output_var} ${${input_var}} PARENT_SCOPE)
  else ()
    set(${output_var} ${ARGN} PARENT_SCOPE)
  endif ()
endfunction ()

# Adds an AMPL solver.
#
# Usage:
#   add_ampl_solver(<name> [DLL_RUNTIME] [MODULE module]
#     [source...] [INCLUDE_DIRS dir...] [LIBRARIES lib...] [EXTERNAL])
#
# If DLL_RUNTIME option is given the solver is linked with DLL MSVC runtime
# on Windows.
#
# The MODULE option specifies the CMake module name used to check for solver
# libraries. If it is omitted, add_ampl_solver looks in support/cmake for
# a module with a name matching Find<name>.cmake (case-insensitive).
#
# Sources are relative to the solver diretory. If no sources are passed to
# add_ampl_solver, it uses the default list of sources <name>.cc, <name>.h.
#
# INCLUDE_DIRS gives a list of solver include directories. If omitted,
# add_ampl_solver uses <module>_INCLUDE_DIRS.
#
# LIBRARIES gives a list of solver libraries. If omitted, add_ampl_solver
# uses <module>_LIBRARIES.
function(add_ampl_solver name)
  cmake_parse_arguments(add_ampl_solver "DLL_RUNTIME;EXTERNAL"
    "MODULE" "INCLUDE_DIRS;LIBRARIES;DEFINITIONS" ${ARGN})

  set(options)
  if (add_ampl_solver_EXTERNAL)
    set(options EXTERNAL)
  endif ()

  check_module(${name} build_solver ${options})
  if (NOT build_solver)
    return ()
  endif ()

  string(TOUPPER ${name} name_upper)

  # Check if the solver is installed.
  if (add_ampl_solver_MODULE)
    set(module ${add_ampl_solver_MODULE})
  else ()
    set(module )
    file(GLOB modules ${MP_CMAKE_MODULE_DIR}/Find*.cmake)
    foreach (m ${modules})
      string(REGEX REPLACE ".*Find(.*).cmake" "\\1" m ${m})
      string(TOUPPER ${m} module_upper)
      if (module_upper STREQUAL ${name_upper})
        set(module ${m})
        break ()
      endif ()
    endforeach ()
  endif ()
  string(TOUPPER "${module}" module_upper)
  if (module)
    find_package(${module})
    if (NOT ${module_upper}_FOUND)
      return ()
    endif ()
  endif ()

  # Handle arguments.
  get_var(source_names add_ampl_solver_UNPARSED_ARGUMENTS ${name}.cc ${name}.h)
  get_var(include_dirs
    add_ampl_solver_INCLUDE_DIRS ${${module_upper}_INCLUDE_DIRS})
  get_var(libraries add_ampl_solver_LIBRARIES)
  if (NOT libraries AND TARGET ${name}-library)
    set(libraries ${name}-library)
  endif ()
  get_var(definitions
    add_ampl_solver_DEFINITIONS ${${module_upper}_DEFINITIONS})
  if (add_ampl_solver_DLL_RUNTIME AND MSVC)
    set(use_dll_runtime True)
  endif ()

  if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${name})
    # Add directory name to sources.
    set(sources )
    foreach (src ${source_names})
      set(sources ${sources} ${name}/${src})
    endforeach ()
    set(dir ${name}/)
  else ()
    set(sources ${source_names})
    set(dir )
  endif ()

  # Add a solver static library to be linked with the shared library and
  # the exectuable.
  set(static_lib ampl${name}-static)
  add_library(${static_lib} STATIC ${sources})
  target_include_directories(${static_lib}
    PUBLIC ${include_dirs} ${PROJECT_SOURCE_DIR}/solvers)
  target_link_libraries(${static_lib} PUBLIC ${libraries})
  target_compile_options(${static_lib} PUBLIC ${definitions})

  if (use_dll_runtime)
    target_link_libraries(${static_lib} PUBLIC asl-dynrt)
  else ()
    target_link_libraries(${static_lib} PUBLIC asl)
  endif ()

  # Add a solver shared library.
  add_library(ampl${name} SHARED ${PROJECT_SOURCE_DIR}/src/solver-c.cc)
  target_link_libraries(ampl${name} ${static_lib})

  # Add a solver executable.
  add_executable(${name} ${dir}main.cc)
  target_link_libraries(${name} ${static_lib})

  install(TARGETS ${name} DESTINATION bin)
endfunction()

check_module(cplex build_cplex)
if (build_cplex)
  add_subdirectory(cplex)
endif ()

set(GECODE_DISABLE_WARNINGS TRUE)
set(GECODE_USE_QT FALSE)
add_ampl_solver(gecode EXTERNAL LIBRARIES gecodeminimodel)
if (TARGET fzn-gecode)
  # Exclude unused Gecode targets from the build.
  set_target_properties(
    fzn-gecode gecodeflatzinc PROPERTIES EXCLUDE_FROM_ALL TRUE)
endif ()

add_ampl_solver(ilogcp DLL_RUNTIME MODULE CPLEX
  concert.cc concert.h ilogcp.cc ilogcp.h ilogcp_date.h
  LIBRARIES cplex-cp ilocplex)
if (TARGET ilogcp)
  target_compile_definitions(amplilogcp-static PRIVATE ILOGCP_NO_VERS)
  target_compile_definitions(ilogcp PRIVATE main1=main)
endif ()

# Copy CPLEX DLL to the target directory.
if (CPLEX_DLL)
  get_filename_component(CPLEX_DLL_NAME "${CPLEX_DLL}" NAME)
  add_custom_command(TARGET ilogcp
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
    "${CPLEX_DLL}" "$<TARGET_FILE_DIR:ilogcp>/${CPLEX_DLL_NAME}")
endif ()

add_subdirectory(jacop)

add_ampl_solver(localsolver)
add_ampl_solver(localsolver-asl
  MODULE LocalSolver LIBRARIES localsolver-library)

# Copy LocalSolver DLL to the target directory.
if (LOCALSOLVER_DLL)
  get_filename_component(LOCALSOLVER_DLL_NAME "${LOCALSOLVER_DLL}" NAME)
  add_custom_command(TARGET localsolver
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
    "${LOCALSOLVER_DLL}"
    "$<TARGET_FILE_DIR:localsolver>/${LOCALSOLVER_DLL_NAME}")
endif ()

add_subdirectory(path)

add_ampl_solver(smpswriter)

add_ampl_solver(ssdsolver)
if (TARGET ssdsolver)
  add_ampl_library(ssd ssdsolver/ssd.cc)
  target_include_directories(ssd PRIVATE .)
endif ()

add_ampl_solver(sulum)
