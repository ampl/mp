var x{1..3} integer >= -2 <= 2;
minimize o: (4)*x[1] + (2)*x[2] + (2)*x[3];
s.t. c1: (x[2] > 1) ==> (x[3] == -1) else (x[2] == -2);
s.t. c2: alldiff(x[2], x[1], x[3]);
