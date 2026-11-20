var x{1..3} integer >= 0 <= 2;
minimize o: (-4)*x[1] + (-2)*x[2] + (-4)*x[3];
s.t. c1: alldiff(x[1], x[2]);
s.t. c2: max(x[3], x[1]) - min(x[2], x[3]) > 1;
