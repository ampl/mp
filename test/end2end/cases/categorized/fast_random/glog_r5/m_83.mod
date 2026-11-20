var x{1..5} integer >= 0 <= 4;
minimize o: (-1)*x[1] + (0)*x[2] + (3)*x[3] + (2)*x[4] + (0)*x[5];
s.t. c1: (x[4] != 0) or (x[3] != 1) or (x[5] != 4);
s.t. c2: alldiff(x[5], x[2], x[3]);
