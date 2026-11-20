var x{1..3} integer >= 0 <= 3;
maximize o: (4)*x[1] + (2)*x[2] + (-4)*x[3];
s.t. c1: max(x[3], x[2]) - min(x[2], x[2]) > 4;
s.t. c2: (x[1] <= 4) ==> (x[3] >= 4);
s.t. c3: max(x[1], x[2]) - min(x[1], x[3]) != 0;
