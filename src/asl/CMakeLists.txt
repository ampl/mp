check_module(asl build_asl)
if (NOT build_asl)
  return ()
endif ()

add_subdirectory(${CMAKE_SOURCE_DIR}/thirdparty/asl
${CMAKE_CURRENT_BINARY_DIR}/asl)


include(CheckSymbolExists)
include(CheckLibraryExists)



# Public ASL headers.
set(ASL_HEADERS aslbuilder.h aslexpr.h aslexpr-visitor.h
  aslproblem.h aslinterface.h)

#add_prefix(ASL_HEADERS solvers/
#  asl.h asl_pfg.h asl_pfgh.h avltree.h funcadd.h getstub.h jacpdim.h nlp.h
#  nlp2.h psinfo.h errchk.h jac2dim.h obj_adj.h)
set(ASL_SOURCES aslbuilder.cc aslexpr.cc aslproblem.cc aslinterface.cc)
check_symbol_exists(mkstemps stdlib.h HAVE_MKSTEMPS)
if (NOT HAVE_MKSTEMPS)
  set(ASL_SOURCES ${ASL_SOURCES} mkstemps.c)
endif ()




add_mp_library(aslmp STATIC ${ASL_SOURCES} ${ASL_HEADERS} 
  LIBRARIES mp asl
  INCLUDE_DIRECTORIES ${ASL1_INCLUDE_DIRS} 
  ${PROJECT_SOURCE_DIR}/src)

if (NOT WIN32)
  target_link_libraries(aslmp m)
endif ()

if (BUILD_SHARED_LIBS)
  add_library(asl-shared $<TARGET_OBJECTS:asl>)
  set_target_properties(asl-shared PROPERTIES
    VERSION ${MP_VERSION} SOVERSION ${MP_VERSION_MAJOR})
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
  target_compile_definitions(aslmp PRIVATE NO_RUSAGE)
endif ()

target_link_libraries(aslmp ${CMAKE_DL_LIBS})

find_package(MATLAB)
if (MATLAB_FOUND AND MATLAB_BUILD)
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
    if (WIN32)
      set(mex_install_dir bin)
    else ()
      set(mex_install_dir lib)
    endif ()
    install(FILES $<TARGET_PROPERTY:${name},FILENAME>
      DESTINATION ${mex_install_dir})
  endforeach ()
endif ()

#add_executable(gjh solvers/examples/gjh.c)
#target_link_libraries(gjh asl)
#add_c_compiler_flags(gjh -Wno-parentheses)

add_ampl_library(simpbit PRIVATE tables/simpbit.c)
add_ampl_library(fullbit PRIVATE tables/fullbit.c)

add_c_compiler_flags(simpbit -Wno-parentheses)
add_c_compiler_flags(fullbit -Wno-parentheses)

add_ampl_library(ampltabl tables/tableproxy.c tables/tableproxyver.h)

find_package(ODBC)
if (ODBC_FOUND)
  add_library(amplodbc STATIC tables/amplodbc.c)
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
message("ASL1_SOURCE_DIR ${ASL1_SOURCE_DIR}")
add_executable(tableproxy
  tables/tableproxy.c ${ASL1_SOURCE_DIR}/printf.c tables/tableproxyver.h)
set_target_properties(tableproxy
  PROPERTIES OUTPUT_NAME tableproxy${TABLEPROXY_SUFFIX})
target_compile_definitions(tableproxy PRIVATE STAND_ALONE=1 QUOTIFY=1)
target_link_libraries(tableproxy asl)
if (MP_CLANG)
  # add_c_compiler_flags cannot detect that GCC doesn't support
  # -Wno-string-plus-int, so add this option for clang only.
  add_c_compiler_flags(tableproxy -Wno-string-plus-int)
endif ()

# Link with system socket libraries.
if (WIN32)
  target_link_libraries(ampltabl wsock32)
  target_link_libraries(tableproxy wsock32)
  target_compile_definitions(tableproxy PRIVATE NO_sockaddr_in6)
elseif (CMAKE_SYSTEM_NAME MATCHES "SunOS")
  target_link_libraries(tableproxy socket nsl)
endif ()

# Suppress warnings in ASL and tableproxy.
if (MSVC)
  foreach (target asl-core asl-core-dynrt gjh ampltabl amplodbc
                  fullbit simpbit tableproxy)
    if (TARGET ${target})
      target_compile_options(${target}
        PRIVATE /wd4013 /wd4018 /wd4101 /wd4244 /wd4273 /wd4267 /wd4996)
    endif ()
  endforeach ()
else ()
  add_c_compiler_flags(aslmp -Wno-unused-result -Wno-parentheses)
endif ()

install(FILES ${ASL_HEADERS} ${CMAKE_CURRENT_BINARY_DIR}/arith.h
              solvers/opcode.hd solvers/r_opn.hd
        DESTINATION include/asl)
install(TARGETS aslmp tableproxy DESTINATION lib RUNTIME DESTINATION bin)
