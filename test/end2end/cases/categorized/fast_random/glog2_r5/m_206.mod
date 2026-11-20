var x{1..6} integer >= 0 <= 2;

minimize o: 0*x[1] + 0*x[2] + 0*x[3] + 0*x[4] + 0*x[5] + 0*x[6];

s.t. c1: (x[1] < 2) ==> (x[1] >= 1);
s.t. c2: ((if x[3] == 2 then 1 else 0) + (if x[3] >= 2 then 1 else 0)) >= 1;
