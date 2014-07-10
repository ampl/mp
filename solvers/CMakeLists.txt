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

if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/amplgsl/gsl/CMakeLists.txt)
  add_subdirectory(amplgsl)
endif ()

add_subdirectory(cp)
add_subdirectory(cplex)

if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/gecode/lib/CMakeLists.txt)
  add_subdirectory(gecode)
endif ()

add_subdirectory(ilogcp)

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

add_subdirectory(path)

if (HAVE_AUTO)
  add_subdirectory(smpswriter)
endif ()

add_subdirectory(ssdsolver)

find_package(Sulum)
if (SULUM_FOUND)
  add_subdirectory(sulum)
endif ()

add_subdirectory(util)

install(TARGETS amplsolver
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  ARCHIVE DESTINATION lib)
