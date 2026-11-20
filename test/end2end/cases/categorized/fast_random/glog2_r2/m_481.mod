var x{1..5} integer >= 0 <= 2;

minimize o: x[1] + 2*x[2] - 5*x[3] - 4*x[4] + 4*x[5];

s.t. c1: atmost 2 {i in 1..5} (x[i] >= -1);
