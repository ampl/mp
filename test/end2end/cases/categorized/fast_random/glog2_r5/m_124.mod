var x{1..5} integer >= -2 <= 3;
minimize o: (4)*x[1] + (-4)*x[2] + (-5)*x[3] + (1)*x[4] + (0)*x[5];
s.t. c1: ((if x[4] != 2 then 1 else 0) + (if x[5] >= 2 then 1 else 0)) >= 1;
s.t. c2: numberof 1 in (x[1], x[2], x[3], x[4], x[5]) == 1;
