var x{1..6} integer >= 0 <= 3;
minimize o: (2)*x[1] + (-2)*x[2] + (-4)*x[3] + (-2)*x[4] + (-3)*x[5] + (6)*x[6];
s.t. c1: (x[4] >= -1) ==> (x[3] >= 3) else (x[4] <= 0);
s.t. c2: abs(x[4] - x[4]) < 1;
s.t. c3: alldiff(x[6], x[5]);
