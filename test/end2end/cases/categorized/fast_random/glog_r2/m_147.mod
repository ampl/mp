var x{1..5} integer >= 0 <= 2;
minimize o: (4)*x[1] + (-3)*x[2] + (4)*x[3] + (-5)*x[4] + (3)*x[5];
s.t. c1: max(x[1], x[1]) - min(x[2], x[3]) <= 1;
