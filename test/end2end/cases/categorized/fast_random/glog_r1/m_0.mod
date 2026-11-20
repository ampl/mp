var x{1..5} integer >= -2 <= 1;
minimize o: (4)*x[1] + (5)*x[2] + (-5)*x[3] + (5)*x[4] + (3)*x[5];
s.t. c1: max(x[2], x[3]) - min(x[3], x[3]) == 1;
s.t. c2: abs(x[5] - x[3]) <= 1;
