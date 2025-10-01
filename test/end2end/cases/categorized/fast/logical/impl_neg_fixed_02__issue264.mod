##############################################
# #264 example 2
##############################################

# Example 2
option presolve 0;

option show_stats 1;
option reset_initial_guesses 1;
option display_precision 0;
option solver gurobi;

param r_p_min = 20;
param r_p_max = 70;

var p binary;
var r;

subj to con1: p = 1 ==> r_p_min <= r <= r_p_max;
subj to con2: p = 1 <== r_p_min <= r <= r_p_max;
subj to fix: p = 0;

