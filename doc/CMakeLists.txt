# Path to the script to extract AMPLGSL docs.
set(EXTRACT_DOCS ${CMAKE_CURRENT_SOURCE_DIR}/extract-docs.py)

add_custom_command(OUTPUT html/index.html
  COMMAND python ${EXTRACT_DOCS}
  COMMAND sphinx-build -b html ${CMAKE_CURRENT_SOURCE_DIR} ../html
  DEPENDS conf.py index.rst ../src/gsl/amplgsl.c ${EXTRACT_DOCS})
add_custom_target(doc DEPENDS html/index.html)
