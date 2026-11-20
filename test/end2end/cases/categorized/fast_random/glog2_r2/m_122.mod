var x{1..6} integer >= -2 <= 2;
minimize o: (-5)*x[1] + (3)*x[2] + (3)*x[3] + (0)*x[4] + (3)*x[5] + (5)*x[6];
s.t. c1: alldiff(x[1], x[5], x[2]);
