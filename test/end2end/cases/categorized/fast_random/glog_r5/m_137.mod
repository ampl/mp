var x{1..6} integer >= -2 <= 1;

minimize o: 2*x[1] + 6*x[2] + 6*x[3] - x[4] + 3*x[5] - 2*x[6];

s.t. c1: ((x[3] >= 0) ==> (x[4] != 0)) and ((x[4] != 0) ==> (x[2] >= 0));
s.t. c2: ((if x[1] > -2 then 1 else 0) + (if x[1] >= -2 then 1 else 0)) >= 1;
