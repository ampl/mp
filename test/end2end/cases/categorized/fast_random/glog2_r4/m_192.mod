var x{1..6} integer >= -2 <= 1;

maximize o: 3*x[1] + 0*x[2] + 5*x[3] + 3*x[4] - 2*x[5] - 2*x[6];

s.t. c1: (x[5] >= 1) ==> (x[4] >= 1);
s.t. c2: ((if x[1] >= 1 then 1 else 0) + (if x[1] >= 1 then 1 else 0)) >= 1;
