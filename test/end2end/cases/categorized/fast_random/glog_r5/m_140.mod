var x{1..4} integer >= -2 <= 2;

minimize o: 5*x[1] - 3*x[2] + 0*x[3] - x[4];

s.t. c1: atmost 1 {i in 1..4} (x[i] >= -2);
