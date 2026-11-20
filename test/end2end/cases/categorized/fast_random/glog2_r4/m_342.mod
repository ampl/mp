var x{1..6} integer >= -2 <= 0;
minimize o: (-5)*x[1] + (5)*x[2] + (0)*x[3] + (-4)*x[4] + (2)*x[5] + (5)*x[6];
s.t. c1: alldiff(x[5], x[6]);
