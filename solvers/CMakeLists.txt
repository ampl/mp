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
#     [source...] [INCLUDE_DIRS dir...] [LIBRARIES lib...])
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
  string(TOUPPER ${name} name_upper)

  cmake_parse_arguments(add_ampl_solver
    DLL_RUNTIME MODULE "INCLUDE_DIRS;LIBRARIES;DEFINITIONS" ${ARGN})

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
  get_var(libraries add_ampl_solver_LIBRARIES ${${module_upper}_LIBRARIES})
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
  if (use_dll_runtime)
    set(sources ${sources} ${PROJECT_SOURCE_DIR}/solvers/util/util-all.cc)
  endif ()

  # Add a solver static library to be linked with the shared library and
  # the exectuable.
  set(static_lib ampl${name}-static)
  add_library(${static_lib} ${sources})
  target_include_directories(${static_lib}
    PUBLIC ${include_dirs} ${PROJECT_SOURCE_DIR}/solvers)
  target_link_libraries(${static_lib} ${libraries})
  target_compile_options(${static_lib} PUBLIC ${definitions})

  if (use_dll_runtime)
    if (NOT TARGET asl-dynrt)
      add_library(asl-dynrt ${ASL_SOURCES})
      target_compile_options(asl-dynrt PUBLIC /MD)
    endif ()
    target_link_libraries(${static_lib} asl-dynrt)
  else ()
    target_link_libraries(${static_lib} asl)
  endif ()

  # Add a solver shared library.
  add_library(ampl${name} SHARED ${PROJECT_SOURCE_DIR}/src/asl/solver-c.cc)
  target_link_libraries(ampl${name} ${static_lib})

  # Add a solver executable.
  add_executable(${name} ${dir}main.cc)
  target_link_libraries(${name} ${static_lib})
endfunction()

add_subdirectory(cplex)

if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/gecode/lib/CMakeLists.txt)
  set(GECODE_DISABLE_WARNINGS True)
  add_subdirectory(gecode/lib)
  add_ampl_solver(gecode INCLUDE_DIRS gecode/lib LIBRARIES gecodeminimodel)
endif ()

find_package(CPLEX)
add_ampl_solver(ilogcp DLL_RUNTIME MODULE CPLEX
  concert.cc concert.h ilogcp.cc ilogcp.h ilogcp_date.h
  INCLUDE_DIRS ${CPLEX_CP_INCLUDE_DIRS} ${CPLEX_ILOCPLEX_INCLUDE_DIRS}
  LIBRARIES ${CPLEX_CP_LIBRARIES} ${CPLEX_ILOCPLEX_LIBRARIES} ${CPLEX_LIBRARIES}
  DEFINITIONS ${CPLEX_CONCERT_DEFINITIONS} ${CPLEX_ILOCPLEX_DEFINITIONS})

# Copy CPLEX DLL to the target directory.
if (CPLEX_DLL)
  get_filename_component(CPLEX_DLL_NAME "${CPLEX_DLL}" NAME)
  add_custom_command(TARGET ilogcp
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
    "${CPLEX_DLL}" "$<TARGET_FILE_DIR:ilogcp>/${CPLEX_DLL_NAME}")
endif ()

add_subdirectory(jacop)

add_ampl_solver(localsolver)
add_subdirectory(path)

if (HAVE_AUTO)
  add_ampl_solver(smpswriter)
endif ()

add_ampl_solver(ssdsolver)
add_ampl_library(ssd ssdsolver/ssd.cc)
target_include_directories(ssd PRIVATE .)

add_ampl_solver(sulum)
