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
  message(STATUS "Target doc disabled (requires sphinx, see requirements.txt)")
  return()
endif()

# A sphinx-build found on PATH may belong to a Python environment without
# the theme and extensions that source/conf.py loads; the doc target then
# fails inside Sphinx ("Could not import extension ...") and, since it is
# part of ALL, the whole build fails. Check them in the interpreter that
# runs sphinx-build (its shebang line), or in the first Python found.
set(MP_DOC_PYTHON_MODULES
  sphinx breathe ampl_sphinx_theme sphinx_design sphinx_sitemap
  sphinxcontrib.googleanalytics)
set(SPHINX_PYTHON "")
file(STRINGS ${SPHINX_EXECUTABLE} _sphinx_shebang LIMIT_COUNT 1)
if (_sphinx_shebang MATCHES "^#!([^ \t]+)")
  set(SPHINX_PYTHON ${CMAKE_MATCH_1})
  if (SPHINX_PYTHON MATCHES "/env$")        # '#!/usr/bin/env python3'
    string(REGEX REPLACE "^#![^ \t]+[ \t]+([^ \t]+).*" "\\1"
           _sphinx_py_name "${_sphinx_shebang}")
    find_program(SPHINX_PYTHON_FROM_ENV ${_sphinx_py_name})
    set(SPHINX_PYTHON ${SPHINX_PYTHON_FROM_ENV})
  endif ()
endif ()
if (NOT SPHINX_PYTHON)
  find_package(Python3 COMPONENTS Interpreter QUIET)
  set(SPHINX_PYTHON ${Python3_EXECUTABLE})
endif ()
if (SPHINX_PYTHON)
  set(_missing_doc_modules "")
  foreach (_mod ${MP_DOC_PYTHON_MODULES})
    execute_process(COMMAND ${SPHINX_PYTHON} -c "import ${_mod}"
                    RESULT_VARIABLE _mod_res OUTPUT_QUIET ERROR_QUIET)
    if (NOT _mod_res EQUAL 0)
      list(APPEND _missing_doc_modules ${_mod})
    endif ()
  endforeach ()
  if (_missing_doc_modules)
    message(STATUS "Target doc disabled: ${SPHINX_PYTHON} (used by "
      "${SPHINX_EXECUTABLE}) lacks Python module(s) ${_missing_doc_modules}. "
      "Install doc/requirements.txt into that environment, or point "
      "SPHINX_EXECUTABLE to a sphinx-build that has them.")
    return()
  endif ()
endif ()

message(STATUS
  "Documentation building configured. "
  "Set BUILD_DOC=off otherwise.")

add_prefix(doc_deps ../
  ${MP_ALL_HEADERS})

set(doc_deps ${doc_deps} ${NLW2_SRC})

set(DOC_BASEDIR ${CMAKE_CURRENT_SOURCE_DIR}/source/)
set(DOC_OUTPUTDIR ${CMAKE_CURRENT_BINARY_DIR}/html/)
set(SPHINX_INDEX_FILE ${DOC_OUTPUTDIR}/index.html)

add_prefix(mp_docs ${DOC_BASEDIR}/
    model-guide.rst
    modeling-overview.rst
    modeling-expressions.rst
    modeling-mo.rst
    modeling-efficiency.rst
    modeling-numerics.rst
    modeling-tools.rst
    modeling-troublesh.rst
    features-guide.rst
    models.rst
    intro.rst components.rst drivers.rst
    developers.rst
    howto.rst
    testing.rst
    details.rst
  )
add_prefix(mp_docs_models ${DOC_BASEDIR}/models/
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
