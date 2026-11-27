var x{1..3} integer >= 0 <= 2;

minimize o: -5*x[1] - 5*x[2] + 2*x[3];

s.t. c1: (x[3] <= 3) and (x[1] != 0) ==> (x[3] >= 3);
s.t. c2: numberof 3 in (x[1], x[2], x[3]) >= 3;
s.t. c3: ((if x[1] >= 0 then 1 else 0) + (if x[3] >= 0 then 1 else 0)) >= 1;
