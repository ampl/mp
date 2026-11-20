var x{1..6} integer >= -2 <= 2;
minimize o: (6)*x[1] + (-2)*x[2] + (3)*x[3] + (-4)*x[4] + (1)*x[5] + (1)*x[6];
s.t. c1: atleast 3 {i in 1..6} (x[i] >= 2);
s.t. c2: max(x[6], x[5]) - min(x[5], x[3]) > 2;
