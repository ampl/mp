var x{1..5} integer >= 0 <= 4;

maximize o: 0*x[1] - 2*x[2] + 0*x[3] + 6*x[4] - 4*x[5];

s.t. c1: ((if x[1] != 1 then 1 else 0) + (if x[5] >= 1 then 1 else 0)) >= 1;
s.t. c2: (x[2] > 3) ==> (x[1] != 2) else (x[5] < 4);
s.t. c3: ((if x[5] < 3 then 1 else 0) + (if x[5] >= 3 then 1 else 0)) >= 1;
