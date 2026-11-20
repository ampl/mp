var x{1..5} integer >= 0 <= 2;

minimize o: 0*x[1] - 3*x[2] - 5*x[3] - x[4] + 4*x[5];

s.t. c1: forall {i in 1..4} x[i] != x[i+1];
s.t. c2: (x[2] >= 0) ==> (x[3] > 0);
