var x{1..6} integer >= 0 <= 2;
minimize o: (2)*x[1] + (1)*x[2] + (3)*x[3] + (3)*x[4] + (-5)*x[5] + (-2)*x[6];
s.t. c1: alldiff(x[1], x[5]);
