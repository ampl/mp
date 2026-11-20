var x{1..3} integer >= -3 <= 1;
maximize o: (-4)*x[1] + (5)*x[2] + (5)*x[3];
s.t. c1: not ((x[1] < -2) and (x[1] < 0));
