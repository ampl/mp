var x{1..3} integer >= -2 <= 1;
maximize o: (-2)*x[1] + (6)*x[2] + (6)*x[3];
s.t. c1: alldiff(x[3], x[1]);
