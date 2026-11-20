var x{1..5} integer >= -2 <= 0;

maximize o: 3*x[1] + 3*x[2] + 0*x[3] - 4*x[4] - 4*x[5];

s.t. c1: (x[3] <= -1) ==> (x[2] > -2);
s.t. c2: ((if x[2] > -1 then 1 else 0) + (if x[2] >= -1 then 1 else 0)) >= 1;
