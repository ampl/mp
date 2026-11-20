var x{1..3} integer >= -2 <= 0;
maximize o: (4)*x[1] + (-5)*x[2] + (4)*x[3];
s.t. c1: max(x[1], x[1]) - min(x[1], x[1]) != -2;
