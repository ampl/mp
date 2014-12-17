# Path to the script to extract AMPLGSL docs.
set(EXTRACT_DOCS ${CMAKE_CURRENT_SOURCE_DIR}/extract-docs.py)

# AMPLGSL documentation source.
add_prefix(amplgsl_docs amplgsl/
  accuracy.rst freedoc.rst history.rst no-warranty.rst examples.rst
  front-matter.rst index.rst rng.rst fdl.rst gpl.rst intro.rst)

# Add commands to copy rst files to the binary directory.
set(doc_depends ${EXTRACT_DOCS})
foreach (doc ${amplgsl_docs} index.rst)
  set(src ${CMAKE_CURRENT_SOURCE_DIR}/${doc})
  set(dst ${CMAKE_CURRENT_BINARY_DIR}/${doc})
  add_custom_command(OUTPUT ${dst}
    COMMAND ${CMAKE_COMMAND} -E copy ${src} ${dst} DEPENDS ${src})
  set(doc_depends ${doc_depends} ${dst})
endforeach ()

add_custom_command(OUTPUT html/index.html
  COMMAND python ${EXTRACT_DOCS}
  COMMAND sphinx-build -b html -c ${CMAKE_CURRENT_SOURCE_DIR} . ../html
  DEPENDS conf.py index.rst ../src/gsl/amplgsl.c ${doc_depends})
add_custom_target(doc DEPENDS html/index.html)
