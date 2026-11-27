var x{1..4} integer >= 0 <= 4;

maximize o: -x[1] + 2*x[2] + 3*x[3] - 3*x[4];

s.t. c1: (x[3] != 1) and (x[4] == 0) ==> (x[1] >= 1);
s.t. c2: numberof 3 in (x[1], x[2], x[3], x[4]) == 3;
s.t. c3: (x[3] >= 4) <==> (x[3] + x[1] != 3);
