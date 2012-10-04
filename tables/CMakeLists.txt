add_library(simpbit SHARED simpbit.c)
set_target_properties(simpbit PROPERTIES PREFIX "")
set_target_properties(simpbit PROPERTIES SUFFIX ".dll")

add_library(fullbit SHARED fullbit.c)
set_target_properties(fullbit PROPERTIES PREFIX "")
set_target_properties(fullbit PROPERTIES SUFFIX ".dll")
