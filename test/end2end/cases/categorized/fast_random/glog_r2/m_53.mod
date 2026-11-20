var x{1..5} integer >= -2 <= 0;
minimize o: (-4)*x[1] + (2)*x[2] + (1)*x[3] + (0)*x[4] + (5)*x[5];
s.t. c1: alldiff(x[2], x[3]);
