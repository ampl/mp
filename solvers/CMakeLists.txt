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
#   add_ampl_solver(<name> [DLL_RUNTIME] [SHARED_LIB] [OLD_SHARED_LIB]
#     [EXENAME exe_name]
#     [MODULE module]
#     [source...] [INCLUDE_DIRS dir...] [LIBRARIES lib...] [EXTERNAL])
#
# If DLL_RUNTIME option is given the solver is linked with DLL MSVC runtime
# on Windows.
#
# If SHARED_LIB option is given then the solver is linked to <name>-lib.c
# to export an AMPLS C interface.
#
# If, alternatively, OLD_SHARED_LIB is given, the solver is
# linked to mp/solver-c.cc which assumes create_<solver> function.
# This is deprecated.
#
# The MODULE option specifies the CMake module name used to check for solver
# libraries. If it is omitted, add_ampl_solver looks in support/cmake for
# a module with a name matching Find<name>.cmake (case-insensitive).
#
# Sources are relative to the solver directory. If no sources are passed to
# add_ampl_solver, it uses the default list of sources <name>.cc, <name>.h.
#
# INCLUDE_DIRS gives a list of solver include directories. If omitted,
# add_ampl_solver uses <module>_INCLUDE_DIRS.
#
# LIBRARIES gives a list of solver libraries. If omitted, add_ampl_solver
# uses <module>_LIBRARIES.
function(add_ampl_solver name)
  cmake_parse_arguments(add_ampl_solver
    "DLL_RUNTIME;EXTERNAL;SHARED_LIB;OLD_SHARED_LIB;OVERRIDELIBS"
    "MODULE;EXENAME;VERSION"
    "ADDITIONAL_SOURCES;INCLUDE_DIRS;LIBRARIES;DEFINITIONS"
    ${ARGN})
  get_property(MP_BASE_DIR GLOBAL PROPERTY MP_BASE_DIR)
  set(options)
  if (add_ampl_solver_EXTERNAL)
    set(options EXTERNAL)
  endif ()

  check_module(${name} build_solver ${options})
  if (NOT build_solver)
    return ()
  endif ()
  
  if(add_ampl_solver_OVERRIDELIBS)
     message("Overriding libs detection, using ${libraries} directly")
  else()
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
      if( (NOT ${module_upper}_INCLUDE_DIRS) AND (NOT ${module_upper}_LIBS))
        find_package(${module})
        if (NOT ${module_upper}_FOUND)
          return ()
        endif ()
      else()
        add_library(${module}-library SHARED  IMPORTED)
        foreach(l ${${module_upper}_LIBS})
        set_property(TARGET ${module}-library PROPERTY
          IMPORTED_LOCATION ${l})
        if (MSVC)
            set_property(TARGET ${module}-library PROPERTY
              IMPORTED_IMPLIB ${l})
        endif()
       target_include_directories(${module}-library INTERFACE ${${module_upper}_INCLUDE_DIRS})
       endforeach()
      endif()
    endif ()
  endif()
  # Handle arguments.
  get_var(source_names add_ampl_solver_UNPARSED_ARGUMENTS ${name}.cc ${name}.h)
  get_var(include_dirs
    add_ampl_solver_INCLUDE_DIRS ${${module_upper}_INCLUDE_DIRS})
  get_var(libraries add_ampl_solver_LIBRARIES)

  if (TARGET ${module}-library)
    get_target_property(TRG_PROP_IMPRT_LOC ${module}-library IMPORTED_LOCATION)
    if (NOT libraries OR TRG_PROP_IMPRT_LOC)
      set(libraries ${module}-library)
    endif()
  endif ()

  get_var(definitions
    add_ampl_solver_DEFINITIONS ${${module_upper}_DEFINITIONS})
  if (add_ampl_solver_DLL_RUNTIME AND MSVC)
    set(use_dll_runtime True)
  endif ()
  if (EXISTS ${MP_BASE_DIR}/solvers/${name})
    # Add directory name to sources.
    set(sources )
    foreach (src ${source_names})
      set(sources ${sources} ${MP_BASE_DIR}/solvers/${name}/${src})
    endforeach ()
    set(dir  ${MP_BASE_DIR}/solvers/${name}/)
  else ()
    set(sources ${source_names})
    set(dir )
  endif ()
  # Add additional sources, if present
  if(add_ampl_solver_ADDITIONAL_SOURCES)
    foreach(s ${add_ampl_solver_ADDITIONAL_SOURCES})
      set(sources ${sources} ${s})
    endforeach() 
  endif()

  # Add a solver static library to be linked with the shared library and
  # the executable.
  set(static_lib ampl${name}-static)
  add_library(${static_lib} STATIC ${sources})
  target_include_directories(${static_lib}
    PUBLIC ${include_dirs} ${MP_BASE_DIR}/solvers)
  target_link_libraries(${static_lib} PUBLIC ${libraries})
  target_compile_options(${static_lib} PUBLIC ${definitions})
  if(add_ampl_solver_VERSION)
      target_compile_definitions(${static_lib} PUBLIC DRIVER_DATE=${add_ampl_solver_VERSION})
  endif()
  set_property(TARGET ${static_lib}  PROPERTY CXX_STANDARD 11)
  if (use_dll_runtime)
    target_link_libraries(${static_lib} PUBLIC mp-dynrt)
  else ()
    target_link_libraries(${static_lib} PUBLIC mp)
  endif ()

  if(${add_ampl_solver_SHARED_LIB})
    set(shared_lib ${name}-lib)
    add_library(${shared_lib} SHARED ${MP_BASE_DIR}/solvers/${name}/${name}-lib.c ${MP_BASE_DIR}/src/solver.cc)
    target_compile_definitions(${shared_lib} PRIVATE MP_SOLVER=${name})
    target_link_libraries(${shared_lib} PRIVATE ${static_lib})
  elseif(${add_ampl_solver_OLD_SHARED_LIB})     ## DEPRECATED
    # Old-style shared library if we are not building a directinterface
    # solver. In which case, skip the shared library as we still are to define
    # a proper interface.
    set(shared_lib ampl${name})
    add_library(${shared_lib} SHARED ${MP_BASE_DIR}/src/solver-c.cc)
    target_compile_definitions(${shared_lib} PRIVATE MP_SOLVER=${name})
    target_link_libraries(${shared_lib} PRIVATE ${static_lib})
  endif()

  # Add a solver executable.
  if(add_ampl_solver_EXENAME)
    set(exename ${add_ampl_solver_EXENAME})
    message("Setting exe " ${add_ampl_solver_EXENAME})
  else()
    set(exename ${name})
    message("Setting exe " ${name})
  endif()
  add_executable(${exename} ${dir}main.cc)
  target_link_libraries(${exename} PRIVATE ${static_lib})

  install(TARGETS ${exename} DESTINATION bin)

  add_to_folder(${name} ampl${name}-static ${shared_lib} ${name})
endfunction()

if (NOT MSVC)
  # Remove -pedantic because it causes warnings about long string literals
  # when compiling the cplex driver.
  string(REPLACE "-pedantic" "" CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
endif ()

check_module(cplex build_cplex)
if (build_cplex)
  add_subdirectory(cplex)
endif ()

# Excludes targets from the default build.
function(exclude_from_build)
  foreach (target ${ARGN})
    if (TARGET ${target})
      set_target_properties(${target} PROPERTIES
        EXCLUDE_FROM_ALL TRUE EXCLUDE_FROM_DEFAULT_BUILD TRUE)
    endif ()
  endforeach ()
endfunction()

set(GECODE_DISABLE_WARNINGS TRUE)
set(GECODE_USE_QT FALSE)
add_ampl_solver(gecode EXTERNAL LIBRARIES
  gecodeminimodel gecodesupport gecodekernel gecodeint gecodesearch)
# Exclude unused Gecode targets from the build.
exclude_from_build(fzn-gecode gecodeflatzinc gecode-test)

add_ampl_solver(ilogcp DLL_RUNTIME MODULE CPLEX
  concert.cc concert.h ilogcp.cc ilogcp.h
  LIBRARIES cplex-cp ilocplex ${CMAKE_DL_LIBS})
if (TARGET ilogcp)
  target_compile_definitions(amplilogcp-static PRIVATE ILOGCP_NO_VERS)
  target_compile_definitions(ilogcp PRIVATE main1=main)
endif ()


## Add a solver assuming the Backend/ModelAPI architecture
function(add_ampl_backend name)
    set(SRC ${name}backend.cc ${name}backend.h
            ${name}common.cc ${name}common.h
            ${name}modelapi.cc ${name}modelapi.h
            ${name}-modelapi-connect.cc
            model-mgr-with-std-pb.cc)
    add_ampl_solver(${name} ${SRC} ${ARGN})
    if (TARGET ${name})
      target_compile_definitions(${name} PRIVATE main1=main)
      if(WIN32)
        target_compile_options(ampl${name}-static PUBLIC /bigobj)
      endif()
    endif()
endfunction()

add_ampl_backend(visitor DLL_RUNTIME)

add_ampl_backend(mosek DLL_RUNTIME MODULE MOSEK 
  LIBRARIES ${MOSEK_LIBS} ${CMAKE_DL_LIBS})
add_ampl_backend(copt DLL_RUNTIME SHARED_LIB MODULE COPT 
  LIBRARIES ${COPT_LIBS} ${CMAKE_DL_LIBS})


add_ampl_backend(cbcmp DLL_RUNTIME SHARED_LIB MODULE CBC 
  LIBRARIES ${CBC_LIBS} ${CMAKE_DL_LIBS})


add_ampl_backend(scip DLL_RUNTIME SHARED_LIB MODULE SCIP 
  LIBRARIES ${SCIP_LIBS} ${CMAKE_DL_LIBS})

add_ampl_backend(highsmp DLL_RUNTIME SHARED_LIB MODULE HIGHS 
  LIBRARIES ${HIGHS_LIBS} ${CMAKE_DL_LIBS})

add_ampl_backend(ortoolsmp DLL_RUNTIME SHARED_LIB MODULE ortoolsmp 
  LIBRARIES ${ortoolsmp_LIBS} ${CMAKE_DL_LIBS})

add_ampl_backend(cplexmp DLL_RUNTIME SHARED_LIB MODULE CPLEX
                      LIBRARIES cplex-library ${CMAKE_DL_LIBS})

add_ampl_backend(xpress DLL_RUNTIME SHARED_LIB MODULE XPRESS
                      LIBRARIES ${xpress_LIBS} ${CMAKE_DL_LIBS})


set(GRBSRC gurobibackend.cc gurobibackend.h
            gurobicommon.cc gurobicommon.h
            gurobimodelapi.cc gurobimodelapi.h
            model-mgr-with-std-pb.cc
            gurobi-modelapi-connect.cc
            gurobi-ampls-c-api.h)

add_ampl_solver(gurobi DLL_RUNTIME SHARED_LIB MODULE GUROBI
  ${GRBSRC}
  LIBRARIES ${GUROBI_LIBRARY} ${CMAKE_DL_LIBS})
if (TARGET gurobi)
  target_compile_definitions(gurobi PRIVATE main1=main)
  if(WIN32)
    target_compile_options(amplgurobi-static PUBLIC /bigobj)
  endif()
endif()

## A simple Gurobi driver to test AMPLS C API
add_ampl_solver(gurobi_ampls DLL_RUNTIME
  gurobi-ampls.c
  INCLUDE_DIRS ${include_dirs} ${MP_BASE_DIR}/solvers
  LIBRARIES amplgurobi-static)


# Copy CPLEX DLL to the target directory.
if (CPLEX_DLL AND (TARGET ilogcp))

  get_filename_component(CPLEX_DLL_NAME "${CPLEX_DLL}" NAME)
  add_custom_command(TARGET ilogcp
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
    "${CPLEX_DLL}" "$<TARGET_FILE_DIR:ilogcp>/${CPLEX_DLL_NAME}")
endif ()

add_subdirectory(jacop)

add_ampl_solver(localsolver)
if (TARGET localsolver)
  target_compile_definitions(
    ampllocalsolver-static PRIVATE MP_DRIVER_DATE=20160218)
endif ()

# Copy LocalSolver DLL to the target directory.
if (LOCALSOLVER_DLL)
  get_filename_component(LOCALSOLVER_DLL_NAME "${LOCALSOLVER_DLL}" NAME)
  add_custom_command(TARGET localsolver
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
    "${LOCALSOLVER_DLL}"
    "$<TARGET_FILE_DIR:localsolver>/${LOCALSOLVER_DLL_NAME}")
endif ()

add_ampl_solver(smpswriter)

if (TARGET ilogcp)
  add_ampl_solver(ssdsolver LIBRARIES amplilogcp-static)
  if (TARGET ssdsolver)
    add_ampl_library(ssd ssdsolver/ssd.cc)
    target_include_directories(ssd PRIVATE . 
    ${CMAKE_SOURCE_DIR}/src/asl/solvers)
    if (HAVE_UNIQUE_PTR)      ## Need because of nested target definition
      target_compile_definitions(ssd PUBLIC MP_USE_UNIQUE_PTR)
    endif()
  endif ()
endif ()



set(ODHSRC ../odh/odhcommon.h)
add_prefix(DIRGRBSRC "../gurobi/" ${GRBSRC})
set(GRBODHSRC gurobiodhbackend.h gurobiodhbackend.cc)
add_ampl_solver(gurobiodh
                ${DIRGRBSRC} ${ODHSRC} ${GRBODHSRC}
                DLL_RUNTIME MODULE GUROBI
                LIBRARIES ${GUROBI_LIBRARY} ${CMAKE_DL_LIBS}
                )
if (TARGET gurobiodh)
  target_compile_definitions(gurobiodh PRIVATE main1=main)
  if(WIN32)
    target_compile_options(amplgurobiodh-static PUBLIC /bigobj)
  endif()
endif()




set(CPLEXSRC cplexmpbackend.cc cplexmpbackend.h
            cplexmpcommon.cc cplexmpcommon.h
            cplexmpmodelapi.cc cplexmpmodelapi.h
            cplexmp-modelapi-connect.cc
            model-mgr-with-std-pb.cc)
add_prefix(DIRCPLEXSRC "../cplexmp/" ${CPLEXSRC})
set(CPLEXODHSRC cplexodhbackend.h cplexodhbackend.cc)
add_ampl_solver(cplexodh
                ${DIRGRBSRC} ${CPLEXODHSRC} ${ODHSRC}
                DLL_RUNTIME MODULE CPLEX
                LIBRARIES ${CPLEX_LIBRARY} ${CMAKE_DL_LIBS}
                )
if (TARGET cplexodh)
  target_compile_definitions(cplexodh PRIVATE main1=main)
  if(WIN32)
    target_compile_options(amplcplexodh-static PUBLIC /bigobj)
  endif()
endif()

  
