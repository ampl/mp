var x{1..3} integer >= 0 <= 4;

maximize o: 17*x[1] + 3000*x[2] - 2000*x[3];

s.t. c1: ((if x[3] >= 3 then 1 else 0) + (if x[2] >= 3 then 1 else 0)) >= 1;
s.t. c2: (x[2] >= 4) ==> (x[1] <= 1) else (x[3] <= 0);
