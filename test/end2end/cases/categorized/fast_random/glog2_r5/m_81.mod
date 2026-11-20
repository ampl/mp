var x{1..6} integer >= 0 <= 2;
maximize o: (-1)*x[1] + (-4)*x[2] + (6)*x[3] + (-2)*x[4] + (-3)*x[5] + (-4)*x[6];
s.t. c1: ((if x[4] >= -1 then 1 else 0) + (if x[6] >= -1 then 1 else 0)) >= 1;
s.t. c2: atleast 3 {i in 1..6} (x[i] >= 2);
s.t. c3: numberof 0 in (x[1], x[2], x[3], x[4], x[5], x[6]) > 3;
