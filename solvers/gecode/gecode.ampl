# Declarations of suffixes and parameters for Gecode.

# Constraint suffix that specifies consistency level for integer propagators.
suffix icl integer >= 0 <= 3 IN;

# Possible values for the icl suffix.
param icl_val = 0;
param icl_bnd = 1;
param icl_dom = 2;
param icl_def = 3;
