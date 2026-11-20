var x{1..6} integer >= 0 <= 3;

minimize o: -2*x[1] + 4*x[2] + x[3] - 4*x[4] + 0*x[5] - 5*x[6];

s.t. c1: forall {i in 1..5} x[i] != x[i+1];
s.t. c2: atmost 2 {i in 1..6} (x[i] >= 2);
s.t. c3: not ((x[5] < 1) and (x[1] == 3));
