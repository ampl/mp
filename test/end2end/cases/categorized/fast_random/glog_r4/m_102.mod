var x{1..4} integer >= 0 <= 2;
minimize o: (5)*x[1] + (4)*x[2] + (0)*x[3] + (5)*x[4];
s.t. c1: ((if x[2] >= 2 then 1 else 0) + (if x[3] >= 2 then 1 else 0)) >= 1;
s.t. c2: (x[4] > 2) ==> (x[2] > 0);
s.t. c3: numberof 1 in (x[1], x[2], x[3], x[4]) >= 1;
