var x{1..3} integer >= 0 <= 3;
maximize o: (-2)*x[1] + (0)*x[2] + (-4)*x[3];
s.t. c1: (x[2] == 0) or (x[3] > 0) or (x[1] > 2);
s.t. c2: alldiff(x[3], x[1], x[2]);
s.t. c3: (if x[3] < 1 then x[3] else x[3]) <= 0;
