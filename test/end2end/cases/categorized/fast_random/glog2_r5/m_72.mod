var x{1..6} integer >= -2 <= 2;
minimize o: (1)*x[1] + (-5)*x[2] + (1)*x[3] + (6)*x[4] + (5)*x[5] + (5)*x[6];
s.t. c1: ((if x[6] > 1 then 1 else 0) + (if x[6] >= 1 then 1 else 0)) >= 1;
s.t. c2: numberof 0 in (x[1], x[2], x[3], x[4], x[5], x[6]) < 1;
s.t. c3: (x[6] <= 1) ==> (x[6] < 1) else (x[6] <= 1);
