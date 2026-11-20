var x{1..3} integer >= 0 <= 4;

minimize o: 3*x[1] + 3*x[2] - 3*x[3];

s.t. c1: (x[1] <= 4) or (x[3] == 0) or (x[2] <= 1);
s.t. c2: atmost 2 {i in 1..3} (x[i] >= 0);
