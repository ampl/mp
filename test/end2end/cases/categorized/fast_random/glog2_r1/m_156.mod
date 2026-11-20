var x{1..4} integer >= -2 <= 0;
minimize o: (0)*x[1] + (0)*x[2] + (0)*x[3] + (0)*x[4];
s.t. c1: alldiff(x[3], x[2], x[1], x[4]);
