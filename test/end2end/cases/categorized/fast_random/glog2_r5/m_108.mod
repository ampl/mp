var x{1..4} integer >= 0 <= 1;
minimize o: (4)*x[1] + (6)*x[2] + (0)*x[3] + (-4)*x[4];
s.t. c1: (x[3] != 1) and (x[2] != 0) ==> (x[3] >= 1);
s.t. c2: numberof 1 in (x[1], x[2], x[3], x[4]) == 1;
s.t. c3: ((if x[4] != 1 then 1 else 0) + (if x[1] >= 1 then 1 else 0)) >= 1;
