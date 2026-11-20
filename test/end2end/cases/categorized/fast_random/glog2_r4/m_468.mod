var x{1..4} integer >= 0 <= 1;

maximize o: -3*x[1] + 5*x[2] - 4*x[3] + 0*x[4];

s.t. c1: abs(x[2] - x[4]) < 1;
s.t. c2: (x[4] > 1) and (x[4] == 0) ==> (x[2] >= 1);
s.t. c3: abs(x[2] - x[2]) > 1;
s.t. c4: ((if x[2] > 1 then 1 else 0) + (if x[1] >= 1 then 1 else 0)) >= 1;
