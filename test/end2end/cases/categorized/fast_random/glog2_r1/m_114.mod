var x{1..3} integer >= -3 <= 2;
maximize o: (-2)*x[1] + (5)*x[2] + (0)*x[3];
s.t. c1: alldiff(x[3], x[2], x[1]);
s.t. c2: abs(x[3] - x[3]) > 3;
s.t. c3: not ((x[2] < -2) and (x[3] == -2));
