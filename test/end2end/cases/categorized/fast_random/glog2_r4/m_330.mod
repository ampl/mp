var x{1..6} integer >= 0 <= 2;

minimize o: -3*x[1] - 2*x[2] + 2*x[3] + 0*x[4] - 4*x[5] - 4*x[6];

s.t. c1: atmost 2 {i in 1..6} (x[i] >= 0);
s.t. c2: (x[2] == 2) <==> (x[4] + x[4] < 2);
