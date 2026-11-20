var x{1..3} integer >= 0 <= 3;
minimize o: (1)*x[1] + (0)*x[2] + (0)*x[3];
s.t. c1: alldiff(x[3], x[2], x[1]);
s.t. c2: (x[3] >= 3) ==> (x[2] == 2);
