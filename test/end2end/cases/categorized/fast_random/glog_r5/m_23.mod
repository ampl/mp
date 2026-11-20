var x{1..3} integer >= 0 <= 4;

minimize o: -3*x[1] + 3*x[2] - x[3];

s.t. c1: atmost 2 {i in 1..3} (x[i] >= 0);
