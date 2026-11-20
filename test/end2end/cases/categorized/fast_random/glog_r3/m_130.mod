var x{1..6} integer >= 0 <= 2;
minimize o: (-5)*x[1] + (3)*x[2] + (0)*x[3] + (-5)*x[4] + (6)*x[5] + (1)*x[6];
s.t. c1: (x[2] <= 1) <==> (x[2] + x[1] >= 1);
s.t. c2: max(x[3], x[6]) - min(x[2], x[6]) < 0;
s.t. c3: numberof 1 in (x[1], x[2], x[3], x[4], x[5], x[6]) != 2;
s.t. c4: exists {i in 1..6} x[i] == 2;
