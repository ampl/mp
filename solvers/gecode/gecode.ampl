# Declarations of suffixes and parameters for Gecode.

# Constraint suffix that specifies propagation level for integer propagators.
suffix ipl integer >= 0 <= 3 IN;

# Possible values for the ipl suffix.
param ipl_def = 0;
param ipl_val = 1;
param ipl_bnd = 2;
param ipl_dom = 3;
