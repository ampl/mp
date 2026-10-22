# Per-objective time limit and MIP gap, as accepted by
# every MIP solver and by MP2NL (which passes them on
# to the subsolver under its own keywords).
suffix option_timelimit;
suffix option_mipgap;
suffix objpriority;

let total_cost["A&P"].objpriority := 2;
let total_cost["JEWEL"].objpriority := 1;
let total_cost["VONS"].objpriority := 1;
let total_number.objpriority := 1;

let total_cost["A&P"].option_timelimit := 12;
let total_cost["VONS"].option_timelimit := 3.0;
let total_number.option_timelimit := 3;

let total_cost["A&P"].option_mipgap := 0.0001;
let total_cost["JEWEL"].option_mipgap := 0.001;
