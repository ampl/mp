var x{1..3} integer >= -3 <= 2;
minimize o: (-5)*x[1] + (3)*x[2] + (3)*x[3];
s.t. c1: (x[2] != 0) ==> (x[2] <= -2) else (x[2] >= -3);
s.t. c2: alldiff(x[1], x[2], x[3]);
s.t. c3: not ((x[3] < 1) and (x[1] < 0));
