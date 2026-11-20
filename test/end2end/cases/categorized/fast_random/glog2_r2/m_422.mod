var x{1..4} integer >= -2 <= 2;

minimize o: -2*x[1] - 5*x[2] + 0*x[3] + 4*x[4];

s.t. c1: (if x[1] < 1 then x[4] else x[1]) != 0;
s.t. c2: numberof -2 in (x[1], x[2], x[3], x[4]) < 3;
s.t. c3: forall {i in 1..3} x[i] != x[i+1];
