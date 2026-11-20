var x{1..5} integer >= 0 <= 3;
maximize o: (-4)*x[1] + (-3)*x[2] + (-5)*x[3] + (-2)*x[4] + (3)*x[5];
s.t. c1: alldiff(x[3], x[2], x[5], x[4]);
