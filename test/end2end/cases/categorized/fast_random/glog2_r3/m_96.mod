var x{1..6} integer >= 0 <= 4;

maximize o: 0*x[1] + 0*x[2] + 0*x[3] + 0*x[4] + 0*x[5] + 0*x[6];

s.t. c1: not ((x[1] >= 2) and (x[6] > 1));
s.t. c2: ((if x[1] != 1 then 1 else 0) + (if x[2] >= 1 then 1 else 0)) >= 1;
s.t. c3: (x[2] >= 0) ==> (x[4] > 2);
