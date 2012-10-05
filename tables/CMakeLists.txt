add_ampl_library(simpbit SHARED simpbit.c)
add_ampl_library(fullbit SHARED fullbit.c)
add_ampl_library(ampltabl SHARED tableproxy.c tableproxyver.h)

if (CMAKE_SIZEOF_VOID_P EQUAL 8)
  set(TABLEPROXY_SUFFIX 64)
else ()
  set(TABLEPROXY_SUFFIX 32)
endif ()
set(TABLEPROXY_NAME tableproxy${TABLEPROXY_SUFFIX})
add_executable(${TABLEPROXY_NAME} tableproxy.c tableproxyver.h)
set_target_properties(${TABLEPROXY_NAME}
  PROPERTIES COMPILE_DEFINITIONS STAND_ALONE=1)
target_link_libraries(${TABLEPROXY_NAME} ${SOLVER_LIBS})
