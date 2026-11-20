var x{1..5} integer >= -2 <= 2;
minimize o: (0)*x[1] + (0)*x[2] + (0)*x[3] + (0)*x[4] + (0)*x[5];
s.t. c1: alldiff(x[2], x[5], x[1], x[4]);
s.t. c2: alldiff(x[3], x[1]);
