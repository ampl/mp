# Simple diet-problem script that gets data from tables.
# This one works with the .tab files written by "diettabw.x".
# With the files written by "dietwrit.x", this script will also
# work with ".tab" changed to ".bit" or -- with the addition of
# a "load simpbit.dll;" command -- to ".abt".

model diet.mod;
table calories IN: [NUTR,FOOD] amt;
table foods IN:
	[FOOD] IN, cost, f_min, f_max;
table nutrients IN:
	[NUTR] IN, n_min, n_max;
read table calories;
read table foods;
read table nutrients;
solve;
display _varname, _var, _var.rc;
