var x{1..5} integer >= -2 <= 2;
maximize o: (-1)*x[1] + (-1)*x[2] + (5)*x[3] + (-5)*x[4] + (1)*x[5];
s.t. c1: not ((x[2] > -2) and (x[3] < 1));
