add_custom_command(OUTPUT html/index.html
  COMMAND python extract-docs.py
  COMMAND sphinx-build -b html . ../html
  DEPENDS conf.py index.rst tables/index.rst tables/mysql.rst
          ../solvers/amplgsl/amplgsl.c extract-docs.py)
add_custom_target(doc DEPENDS html/index.html)
