var x{1..4} integer >= -2 <= 2;

minimize o: -x[1] - 5*x[2] + 4*x[3] + 0*x[4];

s.t. c1: forall {i in 1..3} x[i] != x[i+1];
s.t. c2: (if x[2] == 1 then x[4] else x[4]) >= 2;
s.t. c3: ((if x[4] == -1 then 1 else 0) + (if x[1] >= -1 then 1 else 0)) >= 1;
s.t. c4: max(x[2], x[1]) - min(x[2], x[4]) <= -1;
