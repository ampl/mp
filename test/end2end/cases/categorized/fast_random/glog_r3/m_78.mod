var x{1..3} integer >= -2 <= 2;
minimize o: (0)*x[1] + (6)*x[2] + (-2)*x[3];
s.t. c1: max(x[3], x[2]) - min(x[3], x[3]) > 0;
s.t. c2: max(x[1], x[1]) - min(x[1], x[2]) > 1;
s.t. c3: (x[1] == -1) or (x[2] <= -1) or (x[2] >= 1);
s.t. c4: max(x[3], x[1]) - min(x[1], x[3]) < -2;
