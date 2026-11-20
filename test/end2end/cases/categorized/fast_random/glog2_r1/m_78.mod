var x{1..5} integer >= 0 <= 1;
maximize o: (6)*x[1] + (6)*x[2] + (-5)*x[3] + (-4)*x[4] + (-2)*x[5];
s.t. c1: numberof 1 in (x[1], x[2], x[3], x[4], x[5]) == 1;
s.t. c2: numberof 0 in (x[1], x[2], x[3], x[4], x[5]) > 0;
s.t. c3: max(x[5], x[3]) - min(x[3], x[2]) <= 1;
s.t. c4: (x[5] != 1) ==> (x[2] != 0) else (x[1] == 1);
