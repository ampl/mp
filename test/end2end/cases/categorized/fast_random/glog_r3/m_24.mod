var x{1..3} integer >= -2 <= 1;
minimize o: (-4)*x[1] + (-2)*x[2] + (-4)*x[3];
s.t. c1: alldiff(x[2], x[1]);
s.t. c2: abs(x[3] - x[3]) > 3;
s.t. c3: max(x[3], x[2]) - min(x[3], x[2]) == -2;
