var x{1..4} integer >= 0 <= 4;

maximize o: 5*x[1] + 5*x[2] + 4*x[3] - 3*x[4];

s.t. c1: forall {i in 1..3} x[i] != x[i+1];
s.t. c2: abs(x[4] - x[4]) > 4;
s.t. c3: ((if x[3] > 2 then 1 else 0) + (if x[1] >= 2 then 1 else 0)) >= 1;
