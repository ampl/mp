var x{1..6} integer >= 0 <= 4;

maximize o: 3*x[1] + 6*x[2] - 4*x[3] - 4*x[4] - 3*x[5] + 3*x[6];

s.t. c1: atmost 3 {i in 1..6} (x[i] >= 0);
