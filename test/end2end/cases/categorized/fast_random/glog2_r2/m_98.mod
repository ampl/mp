var x{1..6} integer >= 0 <= 1;

minimize o: -5*x[1] + 0*x[2] - x[3] - 4*x[4] - 5*x[5] + x[6];

s.t. c1: (x[2] >= 0) ==> (x[3] <= 0);
s.t. c2: ((if x[6] == 1 then 1 else 0) + (if x[1] >= 1 then 1 else 0)) >= 1;
