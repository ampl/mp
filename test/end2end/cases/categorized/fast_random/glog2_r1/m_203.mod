var x{1..3} integer >= -2 <= 0;
maximize o: (-3)*x[1] + (-1)*x[2] + (0)*x[3];
s.t. c1: alldiff(x[1], x[3], x[2]);
s.t. c2: not ((x[2] > -2) and (x[3] <= 0));
s.t. c3: abs(x[1] - x[3]) != 0;
s.t. c4: x[3] == 0;
