var x{1..4} integer >= 0 <= 4;

minimize o: x[1] + 3*x[2] + 4*x[3] - 3*x[4];

s.t. c1: forall {i in 1..3} x[i] != x[i+1];
s.t. c2: (x[2] == 0) ==> (x[4] > 5) else (x[2] > 5);
