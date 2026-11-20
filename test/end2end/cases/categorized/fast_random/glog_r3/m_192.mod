var x{1..4} integer >= 0 <= 3;

minimize o: 6*x[1] - 4*x[2] - 5*x[3] + 5*x[4];

s.t. c1: (x[1] != 1) or (x[1] <= 0) or (x[1] >= 3);
s.t. c2: (if x[1] != 2 then x[2] else x[4]) > 3;
s.t. c3: forall {i in 1..3} x[i] != x[i+1];
s.t. c4: max(x[2], x[2]) - min(x[4], x[3]) != 2;
