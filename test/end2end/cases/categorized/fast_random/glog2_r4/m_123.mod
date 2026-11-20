var x{1..6} integer >= 0 <= 3;

minimize o: 2*x[1] + 5*x[2] + 5*x[3] - 3*x[4] + 4*x[5] + 6*x[6];

s.t. c1: forall {i in 1..5} x[i] != x[i+1];
s.t. c2: max(x[5], x[5]) - min(x[3], x[3]) <= 0;
s.t. c3: ((if x[3] <= 3 then 1 else 0) + (if x[5] >= 3 then 1 else 0)) >= 1;
s.t. c4: x[6] == 2;
