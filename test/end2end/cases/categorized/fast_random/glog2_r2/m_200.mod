var x{1..4} integer >= 0 <= 3;

maximize o: -4*x[1] + x[2] + 3*x[3] + 3*x[4];

s.t. c1: numberof 2 in (x[1], x[2], x[3], x[4]) != 1;
s.t. c2: (x[4] != 3) ==> (x[4] < 0) else (x[1] != 0);
s.t. c3: ((if x[1] <= 3 then 1 else 0) + (if x[2] >= 3 then 1 else 0)) >= 1;
