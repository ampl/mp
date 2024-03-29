# CMake build script for the AMPL/JaCoP solver.

check_module(jacop build_jacop EXTERNAL)

if (build_jacop)
  file(GLOB JACOP_JAR_PATH ${PROJECT_SOURCE_DIR}/thirdparty/jacop/jacop-*.jar)
  get_filename_component(jar ${JACOP_JAR_PATH} NAME)
  set(JACOP_JAR ${jar} CACHE STRING "JaCoP JAR file name.")
  if (CMAKE_SIZEOF_VOID_P EQUAL 4)
    # Set CMAKE_SYSTEM_PROCESSOR to make sure FindJNI finds the correct
    # AWT and JVM libraries.
    set(CMAKE_SYSTEM_PROCESSOR i586)
  endif ()
  find_package(JNI)
  find_package(Java)
  if (JNI_FOUND AND JAVA_FOUND)
    if (CMAKE_SYSTEM_NAME STREQUAL "Linux")
      # Add libmawt.so to JNI_LIBRARIES to avoid linkage errors on Ubuntu 10.04.
      find_library(JAVA_MAWT_LIBRARY mawt
        PATHS ${JAVA_AWT_LIBRARY_DIRECTORIES} PATH_SUFFIXES xawt)
      if (JAVA_MAWT_LIBRARY)
        set(JNI_LIBRARIES ${JNI_LIBRARIES} ${JAVA_MAWT_LIBRARY})
      endif ()
    endif ()

    add_ampl_solver(jacop jacop.cc jacop.h java.cc java.h
      INCLUDE_DIRS ${JNI_INCLUDE_DIRS} LIBRARIES ${JNI_LIBRARIES})
    if (MSVC)
      target_link_libraries(jacop PRIVATE -DELAYLOAD:jvm.dll delayimp)
    endif ()
    add_custom_command(TARGET jacop
      COMMAND ${CMAKE_COMMAND} -E copy_if_different
              ${JACOP_JAR_PATH} $<TARGET_FILE_DIR:jacop>)

    # Get jacop version.
    if (JACOP_JAR MATCHES "jacop-(.*)\\.jar")
      set(JACOP_VERSION ${CMAKE_MATCH_1})
    else ()
      message(FATAL_ERROR "Cannot determine JaCoP version.")
    endif ()
    target_compile_definitions(ampljacop-static
      PRIVATE JACOP_VERSION=\"${JACOP_VERSION}\")
    
    include(UseJava)
    set(CMAKE_JAVA_INCLUDE_PATH ${JACOP_JAR_PATH})
    set(CMAKE_JAVA_COMPILE_FLAGS -source 1.7 -target 1.7 -Xlint:-options)
    add_jar(ampljacop-java InterruptingListener.java SolutionListener.java
      OUTPUT_NAME ampljacop OUTPUT_DIR ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
    add_custom_command(TARGET ampljacop-java
      COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/ampljacop.jar
        $<TARGET_FILE_DIR:jacop>)
  endif ()
endif ()
