var x{1..4} integer >= 0 <= 3;

minimize o: 6*x[1] + 5*x[2] - 2*x[3] + 6*x[4];

s.t. c1: not ((x[3] != 0) and (x[2] != 3));
s.t. c2: forall {i in 1..3} x[i] != x[i+1];
s.t. c3: max(x[3], x[1]) - min(x[2], x[2]) >= 3;
s.t. c4: numberof 1 in (x[1], x[2], x[3], x[4]) != 2;
