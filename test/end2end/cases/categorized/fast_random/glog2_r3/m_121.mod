var x{1..3} integer >= 0 <= 2;
maximize o: (4)*x[1] + (2)*x[2] + (-1)*x[3];
s.t. c1: max(x[1], x[2]) - min(x[2], x[3]) >= -1;
s.t. c2: (x[2] != 0) ==> (x[1] != 2) else (x[1] > 3);
s.t. c3: not ((x[3] >= 0) and (x[2] == -1));
