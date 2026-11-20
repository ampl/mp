var x{1..5} integer >= 0 <= 4;
minimize o: (4)*x[1] + (-4)*x[2] + (1)*x[3] + (-5)*x[4] + (-2)*x[5];
s.t. c1: max(x[2], x[1]) - min(x[3], x[1]) == 3;
s.t. c2: max(x[2], x[1]) - min(x[3], x[1]) == 3;
