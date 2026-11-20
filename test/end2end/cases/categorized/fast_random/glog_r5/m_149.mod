var x{1..3} integer >= 0 <= 3;
minimize o: (5)*x[1] + (-2)*x[2] + (5)*x[3];
s.t. c1: alldiff(x[3], x[2]);
s.t. c2: not ((x[2] >= 0) and (x[3] > 3));
s.t. c3: (x[3] >= 1) ==> (x[2] >= 3);
s.t. c4: (x[1] > 1) or (x[3] < 3) or (x[1] != 3);
