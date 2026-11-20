var x{1..6} integer >= 0 <= 4;

minimize o: -x[1] - x[2] - 5*x[3] - 5*x[4] + 0*x[5] + x[6];

s.t. c1: (x[5] != 0) <==> (x[3] + x[3] <= 3);
s.t. c2: ((if x[1] > 1 then 1 else 0) + (if x[2] >= 1 then 1 else 0)) >= 1;
