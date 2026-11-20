var x{1..5} integer >= 0 <= 4;

maximize o: 3*x[1] - x[2] + 4*x[3] + 2*x[4] + 4*x[5];

s.t. c1: forall {i in 1..4} x[i] != x[i+1];
s.t. c2: ((if x[1] != 0 then 1 else 0) + (if x[1] >= 0 then 1 else 0)) >= 1;
