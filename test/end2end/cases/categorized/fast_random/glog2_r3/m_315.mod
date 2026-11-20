var x{1..3} integer >= 0 <= 4;

maximize o: 4*x[1] - 5*x[2] - 5*x[3];

s.t. c1: numberof 4 in (x[1], x[2], x[3]) > 1;
s.t. c2: (x[3] >= 2) ==> (x[2] != 0);
s.t. c3: ((if x[2] != 2 then 1 else 0) + (if x[1] >= 2 then 1 else 0)) >= 1;
s.t. c4: (x[1] != 0) ==> (x[2] == 1);
s.t. c5: (x[2] < 1) and (x[1] > 4) ==> (x[1] >= 1);
s.t. c6: max(x[3], x[2]) - min(x[3], x[1]) >= 0;
