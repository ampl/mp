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
  set_property(TARGET arithchk
    PROPERTY COMPILE_DEFINITIONS NO_FPINIT NO_SSIZE_T)
else ()
  set_target_properties(arithchk PROPERTIES COMPILE_DEFINITIONS ASL_NO_FPINITMT)
  target_link_libraries(arithchk m)
endif ()
add_custom_command(OUTPUT arith.h DEPENDS arithchk COMMAND arithchk > arith.h)

include_directories(solvers ${CMAKE_CURRENT_BINARY_DIR})

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
check_symbol_exists(mkstemps stdlib.h HAVE_MKSTEMPS)
if (HAVE_MKSTEMPS)
  add_definitions(-DMP_HAVE_MKSTEMPS=1)
else ()
  set(ASL_SOURCES ${ASL_SOURCES} mkstemps.c)
endif ()

add_library(asl ${ASL_SOURCES}
  aslbuilder.cc aslbuilder.h expr.cc expr.h problem.cc problem.h
  solver-c.cc solver-c.h solver.cc solver.h)

# Link with librt for clock_gettime (Linux on i386).
find_library(RT_LIBRARY rt)
if (RT_LIBRARY)
  target_link_libraries(asl ${RT_LIBRARY})
endif ()
