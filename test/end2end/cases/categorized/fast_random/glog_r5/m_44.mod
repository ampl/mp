var x{1..3} integer >= 0 <= 3;

maximize o: 0*x[1] - 3*x[2] - 3*x[3];

s.t. c1: ((if x[1] > 0 then 1 else 0) + (if x[1] >= 0 then 1 else 0)) >= 1;
s.t. c2: (x[2] < 2) ==> (x[1] > 3);
