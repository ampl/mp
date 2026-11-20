var x{1..5} integer >= 0 <= 1;
maximize o: (-1)*x[1] + (4)*x[2] + (3)*x[3] + (-3)*x[4] + (1)*x[5];
s.t. c1: alldiff(x[1], x[2]);
s.t. c2: max(x[2], x[5]) - min(x[1], x[4]) >= 0;
s.t. c3: numberof 0 in (x[1], x[2], x[3], x[4], x[5]) >= 2;
s.t. c4: max(x[1], x[5]) - min(x[1], x[1]) < 1;
