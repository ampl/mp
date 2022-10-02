# Check requirements.
# Use the specific version of Sphinx (see requirements.txt)
# which fixes local search and has better C++ support.
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

message(STATUS
  "Documentation building configured. "
  "Set BUILD_DOC=off otherwise.")

add_prefix(doc_deps ../ ${MP_ALL_HEADERS})

set(DOC_BASEDIR ${CMAKE_CURRENT_SOURCE_DIR})
set(DOC_OUTPUTDIR ${CMAKE_CURRENT_BINARY_DIR}/html/)
set(SPHINX_INDEX_FILE ${DOC_OUTPUTDIR}/index.html)

add_prefix(mp_docs ${DOC_BASEDIR}/rst/
    model-guide.rst features-guide.rst
    models.rst
    intro.rst components.rst drivers.rst
    developers.rst
    getting-started.rst
    testing.rst
    details.rst
  )
add_prefix(mp_docs_models ${DOC_BASEDIR}/rst/models/
    dietobj.mod dietobj.dat
    iis.mod infeasible.mod
    nqueens.mod
  )
set(mp_docs ${mp_docs} ${mp_docs_models} ${DOC_BASEDIR}/index.rst)
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
