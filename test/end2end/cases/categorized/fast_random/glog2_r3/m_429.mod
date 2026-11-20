var x{1..3} integer >= 0 <= 4;
maximize o: (6)*x[1] + (0)*x[2] + (-4)*x[3];
s.t. c1: not ((x[1] <= 3) and (x[2] > 3));
s.t. c2: numberof 1 in (x[1], x[2], x[3]) < 3;
s.t. c3: (x[3] == 1) and (x[2] > 1) ==> (x[3] >= 1);
s.t. c4: not ((x[1] <= 3) and (x[2] > 3));
s.t. c5: (x[3] == 1) and (x[2] > 1) ==> (x[3] >= 1);
