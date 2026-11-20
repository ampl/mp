var x{1..4} integer >= 0 <= 2;
maximize o: (-4)*x[1] + (-2)*x[2] + (-5)*x[3] + (4)*x[4];
s.t. c1: (x[2] > 1) and (x[2] == 2) ==> (x[2] >= 1);
s.t. c2: numberof 1 in (x[1], x[2], x[3], x[4]) < 2;
s.t. c3: (if x[3] > 2 then x[2] else x[1]) >= 2;
