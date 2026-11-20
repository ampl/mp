var x{1..5} integer >= -2 <= 1;
minimize o: (4)*x[1] + (-4)*x[2] + (2)*x[3] + (2)*x[4] + (2)*x[5];
s.t. c1: max(x[4], x[4]) - min(x[4], x[4]) >= 0;
