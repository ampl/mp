var x{1..5} integer >= 0 <= 2;

maximize o: -5*x[1] - 4*x[2] + 6*x[3] + 3*x[4] + x[5];

s.t. c1: ((if x[4] <= 1 then 1 else 0) + (if x[5] >= 1 then 1 else 0)) >= 1;
s.t. c2: (x[2] >= 2) and (x[3] != 1) ==> (x[3] >= 2);
s.t. c3: ((if x[4] <= 1 then 1 else 0) + (if x[5] >= 1 then 1 else 0)) >= 1;
s.t. c4: (x[2] >= 2) and (x[3] != 1) ==> (x[3] >= 2);
