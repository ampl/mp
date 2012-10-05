add_ampl_library(simpbit simpbit.c)
add_ampl_library(fullbit fullbit.c)

find_package(ODBC)
if (ODBC_FOUND)
  add_library(amplodbc amplodbc.c)
  set(AMPLODBC amplodbc ${ODBC_LIBRARIES})
  set_target_properties(amplodbc
    PROPERTIES COMPILE_DEFINITIONS funcadd_ASL=funcaddodbc_ASL)
  if (HAVE_FPIC_FLAG)
    set_target_properties(amplodbc PROPERTIES COMPILE_FLAGS -fPIC)
  endif ()
endif ()

if (WIN32)
  set(TABLEPROXY_LIBS wsock32)
  add_definitions(-DNO_sockaddr_in6)
endif ()

add_ampl_library(ampltabl tableproxy.c tableproxyver.h)
target_link_libraries(ampltabl ${TABLEPROXY_LIBS} ${AMPLODBC})
set_target_properties(ampltabl
  PROPERTIES COMPILE_DEFINITIONS OTHER_FUNCADD=funcaddodbc_ASL)

if (CMAKE_SIZEOF_VOID_P EQUAL 8)
  set(TABLEPROXY_SUFFIX 64)
else ()
  set(TABLEPROXY_SUFFIX 32)
endif ()
set(TABLEPROXY_NAME tableproxy${TABLEPROXY_SUFFIX})

add_executable(${TABLEPROXY_NAME} tableproxy.c tableproxyver.h)
set_target_properties(${TABLEPROXY_NAME}
  PROPERTIES COMPILE_DEFINITIONS STAND_ALONE=1)
target_link_libraries(${TABLEPROXY_NAME} ${SOLVER_LIBS} ${TABLEPROXY_LIBS})
