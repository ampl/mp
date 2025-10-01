##############################################
# #264 example 1
##############################################

# Example 1
option presolve 10;

option show_stats 1;
option reset_initial_guesses 1;
option display_precision 0;
option solver gurobi;

param r_p_min = 20;
param r_p_max = 70;

var p binary;
var r;

subj to con: p = 1 <==> r_p_min <= r <= r_p_max;
subj to fix: p = 0;

