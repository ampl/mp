var x{1..4} integer >= 0 <= 3;

minimize o: -5*x[1] - 3*x[2] - x[3] - 5*x[4];

s.t. c1: ((if x[3] > 1 then 1 else 0) + (if x[2] >= 1 then 1 else 0)) >= 1;
s.t. c2: forall {i in 1..3} x[i] != x[i+1];
s.t. c3: (x[2] == 2) or (x[2] >= 1) or (x[3] > 3);
