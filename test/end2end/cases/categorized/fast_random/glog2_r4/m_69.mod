var x{1..3} integer >= 0 <= 3;

maximize o: -3*x[1] + 0*x[2] - 3*x[3];

s.t. c1: (x[2] > 2) ==> (x[1] < 2) else (x[1] == 3);
s.t. c2: numberof 0 in (x[1], x[2], x[3]) > 3;
s.t. c3: max(x[1], x[1]) - min(x[1], x[3]) <= 2;
s.t. c4: (x[1] != 1) ==> (x[2] != 3);
s.t. c5: (x[2] > 2) ==> (x[1] < 2) else (x[1] == 3);
s.t. c6: (x[1] != 1) ==> (x[2] != 3);
