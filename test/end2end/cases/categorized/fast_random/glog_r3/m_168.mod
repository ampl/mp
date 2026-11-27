var x{1..5} integer >= -2 <= 0;

minimize o: -3*x[1] - x[2] - 3*x[3] - 5*x[4] - 3*x[5];

s.t. c1: ((x[3] > -1) ==> (x[5] > -2)) and ((x[5] > -2) ==> (x[3] >= -1));
s.t. c2: ((if x[4] == 0 then 1 else 0) + (if x[4] >= 0 then 1 else 0)) >= 1;
