var x{1..4} integer >= 0 <= 3;
maximize o: (-2)*x[1] + (-1)*x[2] + (5)*x[3] + (-3)*x[4];
s.t. c1: (x[1] == 3) and (x[1] <= 3) ==> (x[2] >= 3);
s.t. c2: alldiff(x[4], x[2], x[1], x[3]);
s.t. c3: abs(x[4] - x[2]) == 0;
