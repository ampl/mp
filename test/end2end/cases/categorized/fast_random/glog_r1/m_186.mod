var x{1..5} integer >= 0 <= 3;
maximize o: (3)*x[1] + (1)*x[2] + (5)*x[3] + (0)*x[4] + (-5)*x[5];
s.t. c1: (x[4] >= 1) ==> (x[2] < 1) else (x[4] <= 0);
s.t. c2: not ((x[2] > 1) and (x[2] < 3));
s.t. c3: alldiff(x[2], x[4]);
s.t. c4: numberof 3 in (x[1], x[2], x[3], x[4], x[5]) >= 1;
