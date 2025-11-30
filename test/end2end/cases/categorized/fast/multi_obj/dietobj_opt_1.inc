suffix option_timelimit;
suffix option_mipgap;
suffix option_outlev;      ## int-valued with most solvers
suffix objpriority;

let total_cost["A&P"].objpriority := 2;
let total_cost["JEWEL"].objpriority := 1;
let total_cost["VONS"].objpriority := 1;
let total_number.objpriority := 1;

let total_cost["A&P"].option_timelimit := 12;
let total_cost["JEWEL"].option_timelimit := 0;
let total_cost["VONS"].option_timelimit := 3.0;
let total_number.option_timelimit := 3;
let total_number.option_outlev := 1;


let total_cost["JEWEL"].option_mipgap := 3;


