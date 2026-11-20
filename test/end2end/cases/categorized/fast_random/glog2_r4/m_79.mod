var x{1..3} integer >= -2 <= 2;
maximize o: (-1)*x[1] + (-1)*x[2] + (2)*x[3];
s.t. c1: numberof 0 in (x[1], x[2], x[3]) == 2;
s.t. c2: alldiff(x[1], x[2]);
s.t. c3: numberof 0 in (x[1], x[2], x[3]) == 2;
s.t. c4: alldiff(x[1], x[2]);
