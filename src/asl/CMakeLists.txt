file(READ solvers/details.c0 DETAILS)
string(REPLACE "System_details"
  "${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR}" DETAILS "${DETAILS}")
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/details.c "${DETAILS}")

configure_file(solvers/stdio1.h0 stdio1.h COPYONLY)

# Generate arith.h
if (CMAKE_SYSTEM_PROCESSOR MATCHES "^x86")
  include(CheckTypeSize)
  check_type_size(double DOUBLE_SIZE)
  check_type_size(long LONG_SIZE)
  set(ARITH_H "#define IEEE_8087\n#define Arith_Kind_ASL 1\n")
  set(ARITH_INT ) # Integer type used by arith.h checks.
  math(EXPR LONG_SIZE_X2 "${LONG_SIZE} * 2")
  if (DOUBLE_SIZE EQUAL LONG_SIZE_X2)
    set(ARITH_INT long)
    set(ARITH_INT_SIZE ${LONG_SIZE})
  else ()
    check_type_size(int INT_SIZE)
    math(EXPR INT_SIZE_X2 "${INT_SIZE} * 2")
    if (DOUBLE_SIZE EQUAL INT_SIZE_X2)
      set(ARITH_H "${ARITH_H}#define Long int\n#define Intcast (int)(long)\n")
      set(ARITH_INT int)
      set(ARITH_INT_SIZE ${INT_SIZE})
    endif ()
  endif ()
  if (ARITH_INT)
    check_type_size("struct { double d\; ${ARITH_INT} L\; }[2]" STRUCT_SIZE)
    math(EXPR DOUBLE_PLUS_INT_SIZE_X2 "2 * (${DOUBLE_SIZE} + ${ARITH_INT_SIZE})")
    if (STRUCT_SIZE GREATER DOUBLE_PLUS_INT_SIZE_X2)
      set(ARITH_H "${ARITH_H}#define Double_Align\n")
    endif ()
  endif ()
  if (CMAKE_SIZEOF_VOID_P EQUAL 8)
    set(ARITH_H "${ARITH_H}#define X64_bit_pointers\n")
  endif ()
  # Check long long.
  check_type_size("long long" LONG_LONG_SIZE)
  if (LONG_LONG_SIZE GREATER LONG_SIZE AND
      LONG_LONG_SIZE EQUAL CMAKE_SIZEOF_VOID_P)
    set(ARITH_H "${ARITH_H}#define LONG_LONG_POINTERS\n")
  endif ()
  if (LONG_LONG_SIZE LESS 8)
    set(ARITH_H "${ARITH_H}#define NO_LONG_LONG\n")
  endif ()
  # Check ssize_t.
  check_type_size(size_t SIZE_T_SIZE)
  check_type_size(ssize_t SSIZE_T_SIZE)
  if (NOT SSIZE_T_SIZE)
    if (SIZE_T_SIZE EQUAL LONG_SIZE)
      set(ARITH_SSIZE_T long)
    elseif (SIZE_T_SIZE EQUAL INT_SIZE)
      set(ARITH_SSIZE_T int)
    elseif (SIZE_T_SIZE EQUAL LONG_LONG_SIZE)
      set(ARITH_SSIZE_T "long long")
    else ()
      set(ARITH_SSIZE_T "signed size_t")
    endif ()
    set(ARITH_H "${ARITH_H}#define ssize_t ${ARITH_SSIZE_T}\n")
  elseif (NOT SIZE_T_SIZE EQUAL SSIZE_T_SIZE)
    set(ARITH_H "${ARITH_H}/* sizeof(size_t) = ${SIZE_T_SIZE}")
    set(ARITH_H "${ARITH_H} but sizeof(ssize_t) = ${SSIZE_T_SIZE} */\n")
  endif ()
  file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/arith.h
       "${ARITH_H}#define QNaN0 0x0\n#define QNaN1 0xfff80000\n")
else ()
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
  add_custom_command(OUTPUT arith.h DEPENDS arithchk
    COMMAND ${WINE} $<TARGET_FILE:arithchk> > arith.h)
  if (NOT MSVC)
    string(REPLACE "-Wall -Wextra -pedantic" ""
      CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
    target_compile_options(arithchk PUBLIC -Wno-format-security)
  endif ()
endif ()

set(ASL_CORE_SOURCES
  ${CMAKE_CURRENT_BINARY_DIR}/details.c ${CMAKE_CURRENT_BINARY_DIR}/arith.h)
add_prefix(ASL_CORE_SOURCES solvers/
  errchk.h
  jac2dim.h
  jacpdim.h
  obj_adj.h
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
  set(ASL_CORE_SOURCES ${ASL_CORE_SOURCES} mkstemps.c)
endif ()

if (CMAKE_SIZEOF_VOID_P EQUAL 8)
  # Changing the floating point precision is not supported on x64.
  set(ASL_COMPILE_DEFINITIONS No_Control87)
endif ()

check_symbol_exists(strtoull stdlib.h HAVE_STRTOULL)
if (NOT HAVE_STRTOULL)
  check_symbol_exists(_strtoui64 stdlib.h HAVE_STRTOUI64)
  if (HAVE_STRTOUI64)
    set(ASL_COMPILE_DEFINITIONS ${ASL_COMPILE_DEFINITIONS} strtoull=_strtoui64)
    set(MP_MEX_OPTIONS "-Dstrtoull=_strtoui64")
  endif ()
endif ()

# Compile C sources separately to be able to specify compile options.
add_mp_library(asl-core OBJECT ${ASL_CORE_SOURCES}
  COMPILE_DEFINITIONS ${ASL_COMPILE_DEFINITIONS}
  INCLUDE_DIRECTORIES ${CMAKE_CURRENT_BINARY_DIR})

# Suppress warnings.
if (MSVC)
  foreach (suffix "" "-dynrt")
    target_compile_options(asl-core${suffix} PRIVATE /wd4244 /wd4267)
  endforeach ()
else ()
  add_c_compiler_flags(asl-core -Wno-unused-result -Wno-parentheses)
endif ()

# Public ASL headers.
set(ASL_HEADERS aslbuilder.h aslsolver.h aslexpr.h aslexpr-visitor.h
  aslproblem.h ${CMAKE_CURRENT_BINARY_DIR}/arith.h
${CMAKE_CURRENT_BINARY_DIR}/stdio1.h)
add_prefix(ASL_HEADERS solvers/
  asl.h asl_pfg.h asl_pfgh.h avltree.h funcadd.h getstub.h nlp.h nlp2.h psinfo.h
  opcode.hd r_opn.hd)
set(ASL_SOURCES aslbuilder.cc aslsolver.cc aslexpr.cc aslproblem.cc)

add_mp_library(asl
  ${ASL_SOURCES} ${ASL_HEADERS} OBJECT_LIBRARIES asl-core LIBRARIES mp
  INCLUDE_DIRECTORIES solvers ${PROJECT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_BINARY_DIR})

if (NOT WIN32)
  target_link_libraries(asl m)
endif ()
check_library_exists(dl dlopen "" HAVE_LIBDL)
if (HAVE_LIBDL)
  target_link_libraries(asl dl)
endif ()

if (BUILD_SHARED_LIBS)
  add_library(asl-shared $<TARGET_OBJECTS:asl-core>)
  set_target_properties(asl-shared PROPERTIES OUTPUT_NAME asl)
  target_link_libraries(asl-shared asl)
  install(TARGETS asl-shared DESTINATION lib RUNTIME DESTINATION bin)
endif ()

check_symbol_exists(mkstemps stdlib.h HAVE_MKSTEMPS)
if (HAVE_MKSTEMPS)
  target_compile_definitions(asl PRIVATE MP_HAVE_MKSTEMPS=1)
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

find_package(MATLAB)
if (MATLAB_FOUND)
  set(matlab_asl asl)
  if (MSVC)
    set(matlab_asl asl-dynrt)
  endif ()
  set(solvers_dir ${CMAKE_CURRENT_SOURCE_DIR}/solvers)
  foreach (name amplfunc spamfunc)
    add_mex(${name}
      ${solvers_dir}/examples/${name}.c ${solvers_dir}/funcadd.c
      COMPILE_FLAGS -I${CMAKE_CURRENT_BINARY_DIR}
                    -I${solvers_dir} ${MP_MEX_OPTIONS}
      LIBRARIES ${matlab_asl})
  endforeach ()
endif ()

add_executable(gjh solvers/examples/gjh.c)
target_link_libraries(gjh asl)
add_c_compiler_flags(gjh -Wno-parentheses)

add_ampl_library(simpbit PRIVATE tables/simpbit.c)
add_ampl_library(fullbit PRIVATE tables/fullbit.c)

add_c_compiler_flags(simpbit -Wno-parentheses)
add_c_compiler_flags(fullbit -Wno-parentheses)

add_ampl_library(ampltabl tables/tableproxy.c tables/tableproxyver.h)

find_package(ODBC)
if (ODBC_FOUND)
  add_library(amplodbc STATIC
    ${CMAKE_CURRENT_BINARY_DIR}/arith.h tables/amplodbc.c)
  target_include_directories(amplodbc PRIVATE ${ODBC_INCLUDE_DIR})
  target_compile_definitions(amplodbc PRIVATE funcadd_ASL=funcaddodbc_ASL)
  add_c_compiler_flags(amplodbc -Wno-deprecated-declarations)
  target_link_libraries(amplodbc asl ${ODBC_LIBRARIES})
  target_link_libraries(ampltabl amplodbc)
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

install(FILES ${ASL_HEADERS} DESTINATION include/asl)
install(TARGETS asl tableproxy DESTINATION lib RUNTIME DESTINATION bin)
