var x{1..3} integer >= -2 <= 3;
maximize o: (4)*x[1] + (2)*x[2] + (-5)*x[3];
s.t. c1: not ((x[2] == 0) and (x[1] != -2));
