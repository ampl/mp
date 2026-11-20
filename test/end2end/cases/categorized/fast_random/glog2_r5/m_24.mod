var x{1..4} integer >= -2 <= 2;

minimize o: -x[1] + 0*x[2] + 0*x[3] - 2*x[4];

s.t. c1: ((x[1] <= 1) ==> (x[3] == 1)) and ((x[3] == 1) ==> (x[1] >= 1));
s.t. c2: ((if x[3] >= 0 then 1 else 0) + (if x[3] >= 0 then 1 else 0)) >= 1;
