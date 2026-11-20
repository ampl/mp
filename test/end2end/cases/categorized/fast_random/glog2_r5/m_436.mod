var x{1..6} integer >= -2 <= 0;

minimize o: 0*x[1] + 0*x[2] + 0*x[3] + 0*x[4] + 0*x[5] + 0*x[6];

s.t. c1: forall {i in 1..5} x[i] != x[i+1];
s.t. c2: (x[6] <= -2) ==> (x[5] == 0);
