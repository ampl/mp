var x{1..6} integer >= 0 <= 4;

minimize o: 0*x[1] + 0*x[2] + 0*x[3] + 0*x[4] + 0*x[5] + 0*x[6];

s.t. c1: (x[6] <= 1) and (x[3] <= 3) ==> (x[2] >= 1);
s.t. c2: ((if x[5] > 3 then 1 else 0) + (if x[4] >= 3 then 1 else 0)) >= 1;
s.t. c3: max(x[3], x[6]) - min(x[3], x[6]) <= 0;
