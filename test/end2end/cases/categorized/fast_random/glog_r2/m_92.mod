var x{1..6} integer >= 0 <= 3;
minimize o: (5)*x[1] + (1)*x[2] + (-2)*x[3] + (-3)*x[4] + (6)*x[5] + (5)*x[6];
s.t. c1: max(x[5], x[3]) - min(x[5], x[3]) <= 1;
s.t. c2: abs(x[3] - x[2]) == 2;
s.t. c3: numberof 0 in (x[1], x[2], x[3], x[4], x[5], x[6]) < 2;
s.t. c4: atleast 2 {i in 1..6} (x[i] >= 1);
