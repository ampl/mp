# CMake build script for the AMPL/JaCoP solver.

file(GLOB jar lib/jacop-*.jar)
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

    add_ampl_solver(jacop jacop.cc jacop.h java.cc java.h
      INCLUDE_DIRS ${JNI_INCLUDE_DIRS} LIBRARIES ${JNI_LIBRARIES})
    if (MSVC)
      target_link_libraries(jacop -DELAYLOAD:jvm.dll delayimp)
    endif ()
    add_custom_command(TARGET jacop
      COMMAND ${CMAKE_COMMAND} -E copy_if_different
        lib/${JACOP_JAR} $<TARGET_FILE_DIR:jacop>)

    # Get jacop version.
    if (JACOP_JAR MATCHES "jacop-(.*)\\.jar")
      set(JACOP_VERSION ${CMAKE_MATCH_1})
    else ()
      message(FATAL_ERROR "Cannot determine JaCoP version.")
    endif ()
    target_compile_definitions(ampljacop-static PRIVATE JACOP_VERSION=\"${JACOP_VERSION}\")
    
    include(UseJava)
    set(CMAKE_JAVA_INCLUDE_PATH lib/${JACOP_JAR})
    add_jar(ampljacop-jar Interrupter.java SolutionListener.java
      OUTPUT_NAME ampljacop OUTPUT_DIR ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
    add_custom_command(TARGET ampljacop-jar
      COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/ampljacop.jar $<TARGET_FILE_DIR:jacop>)
  endif ()
endif ()
