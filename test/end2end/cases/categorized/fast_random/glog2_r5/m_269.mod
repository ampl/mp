var x{1..5} integer >= 0 <= 3;
maximize o: (-3)*x[1] + (-1)*x[2] + (-3)*x[3] + (2)*x[4] + (-4)*x[5];
s.t. c1: numberof 0 in (x[1], x[2], x[3], x[4], x[5]) > 1;
s.t. c2: ((if x[4] != 1 then 1 else 0) + (if x[4] >= 1 then 1 else 0)) >= 1;
s.t. c3: numberof 0 in (x[1], x[2], x[3], x[4], x[5]) != 1;
s.t. c4: alldiff(x[5], x[2]);
s.t. c5: numberof 0 in (x[1], x[2], x[3], x[4], x[5]) > 1;
s.t. c6: alldiff(x[5], x[2]);
