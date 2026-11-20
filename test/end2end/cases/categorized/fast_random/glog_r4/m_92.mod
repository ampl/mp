var x{1..6} integer >= 0 <= 3;

minimize o: -5*x[1] + 2*x[2] - 4*x[3] - 2*x[4] + 5*x[5] + 3*x[6];

s.t. c1: atmost 3 {i in 1..6} (x[i] >= 2);
