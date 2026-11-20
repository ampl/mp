var x{1..4} integer >= 0 <= 2;

maximize o: x[1] + 3*x[2] - 4*x[3] + 0*x[4];

s.t. c1: (x[1] != 1) ==> (x[4] > 0) else (x[4] < 0);
s.t. c2: forall {i in 1..3} x[i] != x[i+1];
