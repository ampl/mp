var x{1..6} integer >= 0 <= 1;
minimize o: (2)*x[1] + (0)*x[2] + (5)*x[3] + (1)*x[4] + (-4)*x[5] + (-2)*x[6];
s.t. c1: exists {i in 1..6} x[i] == 0;
s.t. c2: numberof 1 in (x[1], x[2], x[3], x[4], x[5], x[6]) > 0;
s.t. c3: max(x[5], x[5]) - min(x[6], x[5]) >= 1;
