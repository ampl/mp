var x{1..6} integer >= 0 <= 1;

minimize o: -4*x[1] - x[2] + 0*x[3] + 3*x[4] - 2*x[5] - 2*x[6];

s.t. c1: (x[2] >= 0) and (x[3] > 0) ==> (x[3] >= 0);
s.t. c2: atmost 1 {i in 1..6} (x[i] >= 1);
s.t. c3: (x[1] != 1) ==> (x[5] >= 0) else (x[5] == 1);
