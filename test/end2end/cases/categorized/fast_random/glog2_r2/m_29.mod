var x{1..3} integer >= 0 <= 3;
maximize o: (2)*x[1] + (-3)*x[2] + (-4)*x[3];
s.t. c1: alldiff(x[3], x[2]);
s.t. c2: alldiff(x[3], x[2]);
