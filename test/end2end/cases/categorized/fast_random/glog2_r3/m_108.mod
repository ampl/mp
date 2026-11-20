var x{1..4} integer >= 0 <= 1;

minimize o: -5*x[1] + x[2] - 3*x[3] - 5*x[4];

s.t. c1: forall {i in 1..3} x[i] != x[i+1];
s.t. c2: numberof 0 in (x[1], x[2], x[3], x[4]) > 2;
