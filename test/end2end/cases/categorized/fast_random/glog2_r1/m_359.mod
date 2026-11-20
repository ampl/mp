var x{1..4} integer >= 0 <= 4;
maximize o: (5)*x[1] + (-2)*x[2] + (-2)*x[3] + (-5)*x[4];
s.t. c1: alldiff(x[1], x[3], x[2]);
s.t. c2: alldiff(x[2], x[4], x[1]);
s.t. c3: numberof 1 in (x[1], x[2], x[3], x[4]) != 1;
s.t. c4: alldiff(x[1], x[3], x[2]);
s.t. c5: numberof 1 in (x[1], x[2], x[3], x[4]) != 1;
