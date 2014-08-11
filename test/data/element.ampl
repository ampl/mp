load ../solvers/cp/cp.dll;
function element;
var x >= 0 <= 2;
minimize o: element({i in 1..3} i, x);
