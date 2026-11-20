var x{1..4} integer >= 0 <= 3;

maximize o: 6*x[1] - 4*x[2] + 0*x[3] + x[4];

s.t. c1: (x[4] > 2) and (x[4] >= 0) ==> (x[2] >= 2);
s.t. c2: (x[2] >= 0) or (x[3] < 0) or (x[3] == 3);
s.t. c3: (x[1] < 3) <==> (x[2] + x[4] == 1);
s.t. c4: numberof 2 in (x[1], x[2], x[3], x[4]) != 3;
