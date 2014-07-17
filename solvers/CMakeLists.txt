include(CheckFunctionExists)
check_function_exists(getrusage HAVE_GETRUSAGE)
if (NOT HAVE_GETRUSAGE)
 add_definitions(-DNO_RUSAGE)
endif ()

file(READ details.c0 DETAILS)
string(REPLACE "System_details"
  "${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR}" DETAILS "${DETAILS}")
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/details.c "${DETAILS}")

configure_file(stdio1.h0 stdio1.h COPYONLY)

# Generate arith.h
if (NOT WIN32)
  set(FPINIT fpinit.c)
endif ()
add_executable(arithchk arithchk.c ${FPINIT})
if (WIN32)
  set_property(TARGET arithchk
    PROPERTY COMPILE_DEFINITIONS NO_FPINIT NO_SSIZE_T)
else ()
  set_target_properties(arithchk PROPERTIES COMPILE_DEFINITIONS ASL_NO_FPINITMT)
  target_link_libraries(arithchk m)
endif ()
add_custom_command(OUTPUT arith.h DEPENDS arithchk
  COMMAND arithchk > arith.h)
add_custom_target(arith_h SOURCES ${CMAKE_CURRENT_BINARY_DIR}/arith.h)

include(CheckSymbolExists)
check_symbol_exists(strtoull stdlib.h HAVE_STRTOULL)
if (NOT HAVE_STRTOULL)
  check_symbol_exists(_strtoui64 stdlib.h HAVE_STRTOUI64)
  if (HAVE_STRTOUI64)
    add_definitions(-Dstrtoull=_strtoui64)
  endif ()
endif ()

if (CMAKE_SIZEOF_VOID_P EQUAL 8)
  # Changing the floating point precision is not supported on x64.
  add_definitions(-DNo_Control87)
endif ()

set(ASL_SOURCES
  stdio1.h
  asldate.c
  atof.c
  auxinfo.c
  avltree.c
  b_search.c
  basename.c
  bscanf.c
  com2eval.c
  comeval.c
  con1ival.c
  con2ival.c
  con2val.c
  conadj.c
  conpval.c
  conscale.c
  conval.c
  derprop.c
  ${CMAKE_CURRENT_BINARY_DIR}/details.c
  dtoa1.c
  duthes.c
  dynlink.c
  f_read.c
  fg_read.c
  fg_write.c
  fgh_read.c
  fpecatch.c
  fpinit.c
  fullhes.c
  func_add.c
  funcadd1.c
  g_fmt.c
  genrowno.c
  getenv.c
  getstub.c
  htcl.c
  indic_cons.c
  jac0dim.c
  jac2dim.c
  jacdim.c
  jacinc.c
  jacinc1.c
  libnamsave.c
  mach.c
  mainexit.c
  mip_pri.c
  misc.c
  mpec_adj.c
  mypow.c
  names.c
  nl_obj.c
  nqpcheck.c
  obj2val.c
  obj_adj.c
  obj_prec.c
  objconst.c
  objval.c
  objval_.c
  op_type.c
  pfg_read.c
  pfghread.c
  printf.c
  pshvprod.c
  punknown.c
  qp_read.c
  qpcheck.c
  qsortv.c
  readsol.c
  repwhere.c
  rops.c
  rops2.c
  sigcatch.c
  sos_add.c
  sphes.c
  sscanf.c
  stderr.c
  studchk0.c
  suf_sos.c
  value.c
  writesol.c
  wrtsol_.c
  ws_desc.c
  wsu_desc.c
  x2check.c
  xectim.c
  xp1known.c
  xp2known.c)

add_library(amplsolver ${CMAKE_CURRENT_BINARY_DIR}/arith.h ${ASL_SOURCES})
if (NOT WIN32)
  target_link_libraries(amplsolver dl)
endif ()

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
#   add_ampl_solver(<name> [DLL_RUNTIME]
#     [source...] [INCLUDE_DIRS dir...] [LIBRARIES lib...])
#
# Sources are relative to the solver diretory and the default list of sources
# is <name>.cc, <name>.h.
# If DLL_RUNTIME option is given the solver is linked with DLL MSVC runtime
# on Windows.
function(add_ampl_solver name)
  string(TOUPPER ${name} name_upper)

  # Check if the solver is installed.
  set(cmake_dir)
  file(GLOB modules ${AMPL_CMAKE_MODULE_DIR}/Find*.cmake)
  set(module )
  foreach (m ${modules})
    string(REGEX REPLACE ".*Find(.*).cmake" "\\1" m ${m})
    string(TOUPPER ${m} module_upper)
    if (module_upper STREQUAL ${name_upper})
      set(module ${m})
      find_package(${module})
      if (NOT ${module}_FOUND)
        return ()
      endif ()
      break ()
    endif ()
  endforeach ()

  # Parse arguments.
  cmake_parse_arguments(add_ampl_solver
    DLL_RUNTIME "" "INCLUDE_DIRS;LIBRARIES;DEFINITIONS" ${ARGN})
  get_var(source_names add_ampl_solver_UNPARSED_ARGUMENTS ${name}.cc ${name}.h)
  get_var(include_dirs
    add_ampl_solver_INCLUDE_DIRS ${${module}_INCLUDE_DIRS})
  get_var(libraries add_ampl_solver_LIBRARIES ${${module}_LIBRARIES})
  get_var(definitions add_ampl_solver_DEFINITIONS ${${module}_DEFINITIONS})
  if (add_ampl_solver_DLL_RUNTIME AND MSVC)
    set(use_dll_runtime True)
  endif ()

  # Add directory name to sources.
  set(sources )
  foreach (src ${source_names})
    set(sources ${sources} ${name}/${src})
  endforeach ()
  if (use_dll_runtime)
    set(sources ${sources} util/util-all.cc)
  endif ()

  # Add a solver static library to be linked with the shared library and
  # the exectuable.
  set(static_lib ampl${name}-static)
  add_library(${static_lib} ${sources})
  target_include_directories(${static_lib} PUBLIC ${include_dirs})
  target_link_libraries(${static_lib} util ${libraries})
  target_compile_options(${static_lib} PUBLIC ${definitions})
  add_dependencies(${static_lib} arith_h)

  if (use_dll_runtime)
    if (NOT TARGET amplsolver-dynrt)
      add_library(amplsolver-dynrt ${ASL_SOURCES})
      # TODO: set compile flags
      add_dependencies(amplsolver-dynrt arith_h)
    endif ()
    target_link_libraries(${static_lib} amplsolver-dynrt)
  else ()
    target_link_libraries(${static_lib} util)
  endif ()

  # Add a solver shared library.
  add_library(ampl${name} SHARED util/solver-c.cc)
  target_link_libraries(ampl${name} ${static_lib})

  # Add a solver executable.
  add_executable(${name} ${name}/main.cc)
  target_link_libraries(${name} ${static_lib})
endfunction()

if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/amplgsl/gsl/CMakeLists.txt)
  add_subdirectory(amplgsl)
endif ()

add_subdirectory(cp)
add_subdirectory(cplex)

if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/gecode/lib/CMakeLists.txt)
  set(GECODE_DISABLE_WARNINGS True)
  add_subdirectory(gecode/lib)
  add_ampl_solver(gecode INCLUDE_DIRS gecode/lib LIBRARIES gecodeminimodel)
endif ()

find_package(CPLEX)
add_ampl_solver(ilogcp DLL_RUNTIME
  concert.cc concert.h ilogcp.cc ilogcp.h ilogcp_date.h
  INCLUDE_DIRS ${CPLEX_CP_INCLUDE_DIRS} ${CPLEX_ILOCPLEX_INCLUDE_DIRS}
  LIBRARIES ${CPLEX_CP_LIBRARIES} ${CPLEX_ILOCPLEX_LIBRARIES} ${CPLEX_LIBRARIES}
  DEFINITIONS ${CPLEX_CONCERT_DEFINITIONS} ${CPLEX_ILOCPLEX_DEFINITIONS})

# Copy CPLEX DLL to the target directory.
function(copy_cplex_dll target)
  if (CPLEX_DLL)
    get_filename_component(CPLEX_DLL_NAME "${CPLEX_DLL}" NAME)
    add_custom_command(TARGET ${target}
      COMMAND ${CMAKE_COMMAND} -E copy_if_different
      "${CPLEX_DLL}" "$<TARGET_FILE_DIR:${target}>/${CPLEX_DLL_NAME}")
  endif ()
endfunction()
copy_cplex_dll(ilogcp)

file(GLOB jar jacop/lib/jacop-*.jar)
if (jar)
  get_filename_component(jar ${jar} NAME)
  set(JACOP_JAR ${jar} CACHE STRING "JaCoP JAR file name.")
  if (CMAKE_SIZEOF_VOID_P EQUAL 4)
    # Set CMAKE_SYSTEM_PROCESSOR to make sure FindJNI finds the correct
    # AWT and JVM libraries.
    set(CMAKE_SYSTEM_PROCESSOR i586)
  endif ()
  find_package(JNI)
  find_package(Java)
  if (JNI_FOUND AND JAVA_FOUND)
    # Add libmawt.so to JNI_LIBRARIES to avoid linkage errors on Ubuntu 10.04.
    foreach (arch amd64 i386)
      set(LIBMAWT_PATH
        /usr/lib/jvm/java-6-openjdk/jre/lib/${arch}/xawt/libmawt.so)
      if (EXISTS ${LIBMAWT_PATH})
        set(JNI_LIBRARIES ${JNI_LIBRARIES} ${LIBMAWT_PATH})
        break ()
      endif ()
    endforeach ()
    add_subdirectory(jacop)
  endif ()
endif ()

add_ampl_solver(localsolver)
add_subdirectory(path)

if (HAVE_AUTO)
  add_ampl_solver(smpswriter)
endif ()

add_ampl_solver(ssdsolver)
add_ampl_solver(sulum)

add_subdirectory(util)

install(TARGETS amplsolver
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  ARCHIVE DESTINATION lib)
