var x{1..4} integer >= -3 <= 2;

minimize o: -5*x[1] + 0*x[2] - 5*x[3] - x[4];

s.t. c1: (x[1] <= 2) ==> (x[4] != 0) else (x[1] >= 1);
s.t. c2: ((if x[2] >= -2 then 1 else 0) + (if x[2] >= -2 then 1 else 0)) >= 1;
