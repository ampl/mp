var x{1..5} integer >= -3 <= 1;
minimize o: (1)*x[1] + (2)*x[2] + (2)*x[3] + (-4)*x[4] + (6)*x[5];
s.t. c1: alldiff(x[1], x[5], x[2], x[3]);
