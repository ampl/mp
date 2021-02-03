check_module(asl build_asl)
set(build_asl ${build_asl} PARENT_SCOPE)

if (NOT build_asl)
  return()
endif ()
# Check if an asl target is already available (in case mp is used as a submodule
# or we have an imported library). If not, build it.
if(NOT TARGET asl)
  add_subdirectory(${CMAKE_SOURCE_DIR}/thirdparty/asl
    ${CMAKE_CURRENT_BINARY_DIR}/asl)
endif()

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
install(FILES ${ASL_HEADERS} ${CMAKE_CURRENT_BINARY_DIR}/arith.h
              solvers/opcode.hd solvers/r_opn.hd
        DESTINATION include/asl)
install(TARGETS aslmp DESTINATION lib RUNTIME DESTINATION bin)
