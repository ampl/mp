var x{1..3} integer >= 0 <= 3;
maximize o: (5)*x[1] + (-4)*x[2] + (2)*x[3];
s.t. c1: alldiff(x[3], x[1]);
