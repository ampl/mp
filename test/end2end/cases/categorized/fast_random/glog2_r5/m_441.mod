var x{1..3} integer >= 0 <= 4;
maximize o: (-2)*x[1] + (-2)*x[2] + (-5)*x[3];
s.t. c1: numberof 4 in (x[1], x[2], x[3]) > 1;
s.t. c2: (x[3] < 5) and (x[1] <= 4) ==> (x[3] >= 5);
s.t. c3: alldiff(x[1], x[3]);
s.t. c4: (x[1] < 4) and (x[2] != 0) ==> (x[1] >= 4);
