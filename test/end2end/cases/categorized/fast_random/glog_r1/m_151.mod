var x{1..4} integer >= 0 <= 4;
minimize o: (-2)*x[1] + (2)*x[2] + (3)*x[3] + (-1)*x[4];
s.t. c1: abs(x[2] - x[2]) <= 4;
s.t. c2: alldiff(x[2], x[1], x[4], x[3]);
