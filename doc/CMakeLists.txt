# Check requirements.
# Use specific version of Sphinx which fixes local search and has
# better C++ support.
find_program(DOXYGEN doxygen)
if (NOT DOXYGEN)
  message(STATUS "Target doc disabled (requires doxygen)")
  return ()
endif ()
find_package(Sphinx)
if(NOT SPHINX_EXECUTABLE)
  message(STATUS "Documentation generation disabled (requires sphinx)")
  return()
endif()

add_prefix(doc_deps ../ ${MP_ALL_HEADERS})

set(DOC_BASEDIR ${CMAKE_CURRENT_SOURCE_DIR})
set(DOC_OUTPUTDIR ${CMAKE_CURRENT_BINARY_DIR}/html/)
set(SPHINX_INDEX_FILE ${DOC_OUTPUTDIR}/index.html)

add_prefix(mp_docs ${DOC_BASEDIR}/rst/
    common.rst nl-reader.rst flatcvt.rst backend.rst
    reference.rst details.rst model-guide.rst troublesh.rst)
set(mp_docs ${mp_docs} ${DOC_BASEDIR}/index.rst)
add_custom_command(OUTPUT ${SPHINX_INDEX_FILE}
                  COMMAND
                  ${SPHINX_EXECUTABLE} -b html
                  ${DOC_BASEDIR} ${DOC_OUTPUTDIR}
                  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                  DEPENDS ${mp_docs} ${doc_deps}
                  MAIN_DEPENDENCY ${DOC_BASEDIR}/conf.py
                  COMMENT "Generating documentation with Sphinx")

# Named target so we can run the job easily
add_custom_target(mp-doc ALL DEPENDS ${SPHINX_INDEX_FILE}
	SOURCES  ${mp_docs} ${doc_deps})

add_to_folder(docs mp-doc)

install(DIRECTORY ${DOC_OUTPUTDIR}
        DESTINATION docs OPTIONAL
        PATTERN .* EXCLUDE)
