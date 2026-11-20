var x{1..4} integer >= -2 <= 1;

minimize o: -3*x[1] - 3*x[2] + 0*x[3] + 5*x[4];

s.t. c1: atmost 2 {i in 1..4} (x[i] >= 0);
s.t. c2: atmost 2 {i in 1..4} (x[i] >= 0);
