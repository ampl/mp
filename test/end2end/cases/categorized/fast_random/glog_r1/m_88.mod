var x{1..4} integer >= 0 <= 4;
minimize o: (2)*x[1] + (6)*x[2] + (4)*x[3] + (-5)*x[4];
s.t. c1: max(x[2], x[4]) - min(x[3], x[4]) >= 4;
