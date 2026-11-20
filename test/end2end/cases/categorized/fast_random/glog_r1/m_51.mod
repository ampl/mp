var x{1..6} integer >= 0 <= 3;

minimize o: -x[1] + 2*x[2] + 2*x[3] + 5*x[4] + x[5] - 5*x[6];

s.t. c1: forall {i in 1..5} x[i] != x[i+1];
s.t. c2: numberof 3 in (x[1], x[2], x[3], x[4], x[5], x[6]) >= 2;
