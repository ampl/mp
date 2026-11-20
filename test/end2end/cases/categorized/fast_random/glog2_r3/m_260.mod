var x{1..4} integer >= -2 <= 0;
minimize o: (-1)*x[1] + (6)*x[2] + (6)*x[3] + (4)*x[4];
s.t. c1: atmost 2 {i in 1..4} (x[i] >= -1);
