file(READ solvers/details.c0 DETAILS)
string(REPLACE "System_details"
  "${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR}" DETAILS "${DETAILS}")
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/details.c "${DETAILS}")

configure_file(solvers/stdio1.h0 stdio1.h COPYONLY)

# Generate arith.h
if (NOT WIN32)
  set(FPINIT solvers/fpinit.c)
endif ()
add_executable(arithchk solvers/arithchk.c ${FPINIT})
if (WIN32)
  target_compile_definitions(arithchk PRIVATE NO_FPINIT NO_SSIZE_T)
else ()
  target_compile_definitions(arithchk PRIVATE ASL_NO_FPINITMT)
  target_link_libraries(arithchk m)
endif ()
add_custom_command(OUTPUT arith.h DEPENDS arithchk COMMAND arithchk > arith.h)

if (NOT MSVC)
  string(REPLACE "-Wall -Wextra -pedantic" "" CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
  target_compile_options(arithchk PUBLIC -Wno-format-security)
endif ()

set(ASL_SOURCES
  ${CMAKE_CURRENT_BINARY_DIR}/arith.h
  ${CMAKE_CURRENT_BINARY_DIR}/stdio1.h
  ${CMAKE_CURRENT_BINARY_DIR}/details.c)
add_prefix(ASL_SOURCES solvers/
  asl.h
  asl_pfg.h
  asl_pfgh.h
  avltree.h
  errchk.h
  funcadd.h
  getstub.h
  jac2dim.h
  jacpdim.h
  nlp2.h
  nlp.h
  obj_adj.h
  psinfo.h
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

include(CheckSymbolExists)
include(CheckLibraryExists)

check_symbol_exists(mkstemps stdlib.h HAVE_MKSTEMPS)
if (NOT HAVE_MKSTEMPS)
  set(ASL_SOURCES ${ASL_SOURCES} mkstemps.c)
endif ()

add_library(asl-core ${ASL_SOURCES})
target_include_directories(asl-core PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
if (NOT WIN32)
  target_link_libraries(asl-core m)
endif ()
check_library_exists(dl dlopen "" HAVE_LIBDL)
if (HAVE_LIBDL)
  target_link_libraries(asl-core dl)
endif ()

# Suppress warnings.
if (MSVC)
  target_compile_options(asl-core PRIVATE /wd4244 /wd4267)
else ()
  add_c_compiler_flags(asl-core -Wno-unused-result -Wno-parentheses)
endif ()

add_library(asl
  aslbuilder.cc aslbuilder.h aslsolver.cc expr.cc expr.h problem.cc problem.h)
target_include_directories(asl-core PUBLIC solvers ${PROJECT_SOURCE_DIR}/src)
target_link_libraries(asl asl-core mp)

if (CMAKE_SIZEOF_VOID_P EQUAL 8)
  # Changing the floating point precision is not supported on x64.
  target_compile_definitions(asl PRIVATE No_Control87)
endif ()

check_symbol_exists(mkstemps stdlib.h HAVE_MKSTEMPS)
if (HAVE_MKSTEMPS)
  target_compile_definitions(asl PRIVATE MP_HAVE_MKSTEMPS=1)
endif ()

check_symbol_exists(strtoull stdlib.h HAVE_STRTOULL)
if (NOT HAVE_STRTOULL)
  check_symbol_exists(_strtoui64 stdlib.h HAVE_STRTOUI64)
  if (HAVE_STRTOUI64)
    target_compile_definitions(asl-core PRIVATE strtoull=_strtoui64)
  endif ()
endif ()

include(CheckFunctionExists)
check_function_exists(getrusage HAVE_GETRUSAGE)
if (NOT HAVE_GETRUSAGE)
 target_compile_definitions(asl PRIVATE NO_RUSAGE)
endif ()

include(CheckCXXSourceCompiles)
check_cxx_source_compiles("
  #include <unordered_map>
  int main() { std::unordered_map<int, int> m; }" HAVE_UNORDERED_MAP)
if (HAVE_UNORDERED_MAP)
  target_compile_definitions(asl PUBLIC MP_USE_UNORDERED_MAP)
endif ()

find_library(DL_LIBRARY dl)
if (DL_LIBRARY)
  target_link_libraries(asl dl)
endif ()

add_ampl_library(simpbit tables/simpbit.c)
add_ampl_library(fullbit tables/fullbit.c)

add_c_compiler_flags(simpbit -Wno-parentheses)
add_c_compiler_flags(fullbit -Wno-parentheses)

find_package(ODBC)
if (ODBC_FOUND)
  add_library(amplodbc ${CMAKE_CURRENT_BINARY_DIR}/arith.h tables/amplodbc.c)
  target_include_directories(amplodbc PRIVATE ${ODBC_INCLUDE_DIR})
  target_compile_definitions(amplodbc PRIVATE funcadd_ASL=funcaddodbc_ASL)
  add_c_compiler_flags(amplodbc -Wno-deprecated-declarations)
  target_link_libraries(amplodbc asl ${ODBC_LIBRARIES})
  set(AMPLODBC amplodbc ${ODBC_LIBRARIES})
endif ()

add_ampl_library(ampltabl tables/tableproxy.c tables/tableproxyver.h)
target_link_libraries(ampltabl ${AMPLODBC})
if (ODBC_FOUND)
  target_compile_definitions(ampltabl PRIVATE OTHER_FUNCADD=funcaddodbc_ASL)
endif ()

if (CMAKE_SIZEOF_VOID_P EQUAL 8)
  set(TABLEPROXY_SUFFIX 64)
else ()
  set(TABLEPROXY_SUFFIX 32)
endif ()

add_executable(tableproxy
  tables/tableproxy.c solvers/printf.c tables/tableproxyver.h)
set_target_properties(tableproxy
  PROPERTIES OUTPUT_NAME tableproxy${TABLEPROXY_SUFFIX})
target_compile_definitions(tableproxy PRIVATE STAND_ALONE=1 QUOTIFY=1)
target_link_libraries(tableproxy asl)
add_c_compiler_flags(tableproxy -Wno-string-plus-int)

if (WIN32)
  target_link_libraries(ampltabl wsock32)
  target_link_libraries(tableproxy wsock32)
  target_compile_definitions(tableproxy PRIVATE NO_sockaddr_in6)
endif ()
