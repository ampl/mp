var x{1..5} integer >= 0 <= 4;
minimize o: (1)*x[1] + (-1)*x[2] + (4)*x[3] + (3)*x[4] + (1)*x[5];
s.t. c1: max(x[1], x[1]) - min(x[2], x[1]) < 3;
