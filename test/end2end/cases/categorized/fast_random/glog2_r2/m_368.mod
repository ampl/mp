var x{1..5} integer >= 0 <= 1;

maximize o: 0*x[1] - 4*x[2] - 3*x[3] + x[4] - 4*x[5];

s.t. c1: ((if x[2] > 1 then 1 else 0) + (if x[3] >= 1 then 1 else 0)) >= 1;
s.t. c2: (x[1] == 1) and (x[3] >= 0) ==> (x[5] >= 1);
s.t. c3: (x[1] <= 0) <==> (x[2] + x[1] != 0);
