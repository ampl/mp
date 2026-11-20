var x{1..3} integer >= -2 <= 2;

minimize o: -5*x[1] + 6*x[2] - 3*x[3];

s.t. c1: forall {i in 1..2} x[i] != x[i+1];
s.t. c2: (x[2] <= 2) ==> (x[2] >= 1) else (x[2] == 1);
