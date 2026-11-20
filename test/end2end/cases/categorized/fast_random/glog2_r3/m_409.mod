var x{1..5} integer >= 0 <= 3;
minimize o: (1)*x[1] + (5)*x[2] + (0)*x[3] + (1)*x[4] + (3)*x[5];
s.t. c1: alldiff(x[1], x[5], x[2], x[4]);
s.t. c2: alldiff(x[5], x[2], x[4]);
s.t. c3: alldiff(x[1], x[5], x[2], x[4]);
s.t. c4: alldiff(x[5], x[2], x[4]);
