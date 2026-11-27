var x{1..4} integer >= 0 <= 1;

maximize o: -x[1] + 0*x[2] + x[3] - 4*x[4];

s.t. c1: (x[4] != 0) and (x[4] >= 0) ==> (x[2] >= 0);
s.t. c2: numberof 0 in (x[1], x[2], x[3], x[4]) < 3;
s.t. c3: ((if x[2] >= 1 then 1 else 0) + (if x[2] >= 1 then 1 else 0)) >= 1;
