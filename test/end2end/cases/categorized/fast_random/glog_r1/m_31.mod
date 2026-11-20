var x{1..4} integer >= -2 <= 0;

minimize o: -2*x[1] + 4*x[2] - x[3] - 4*x[4];

s.t. c1: atmost 3 {i in 1..4} (x[i] >= -2);
