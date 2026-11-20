var x{1..3} integer >= 0 <= 3;

maximize o: x[1] - 3*x[2] + 5*x[3];

s.t. c1: (x[3] != 2) <==> (x[2] + x[2] == 2);
s.t. c2: (x[3] == 0) and (x[3] <= 0) ==> (x[1] >= 0);
s.t. c3: numberof 2 in (x[1], x[2], x[3]) != 0;
s.t. c4: ((if x[1] >= 0 then 1 else 0) + (if x[1] >= 0 then 1 else 0)) >= 1;
