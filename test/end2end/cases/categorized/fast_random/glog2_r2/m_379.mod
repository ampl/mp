var x{1..3} integer >= 0 <= 2;
minimize o: (6)*x[1] + (3)*x[2] + (-5)*x[3];
s.t. c1: numberof 1 in (x[1], x[2], x[3]) < 2;
s.t. c2: ((if x[3] == 2 then 1 else 0) + (if x[3] >= 2 then 1 else 0)) >= 1;
s.t. c3: numberof 1 in (x[1], x[2], x[3]) < 2;
s.t. c4: ((if x[3] == 2 then 1 else 0) + (if x[3] >= 2 then 1 else 0)) >= 1;
