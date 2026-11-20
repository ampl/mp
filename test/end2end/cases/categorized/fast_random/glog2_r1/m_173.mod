var x{1..6} integer >= 0 <= 3;

maximize o: -3*x[1] + 4*x[2] + 0*x[3] + 4*x[4] - 5*x[5] + x[6];

s.t. c1: (x[4] < 0) and (x[2] != 2) ==> (x[1] >= 0);
s.t. c2: (x[2] > 0) and (x[1] >= 3) ==> (x[1] >= 0);
s.t. c3: numberof 3 in (x[1], x[2], x[3], x[4], x[5], x[6]) <= 2;
s.t. c4: numberof 2 in (x[1], x[2], x[3], x[4], x[5], x[6]) != 0;
s.t. c5: x[6] == 3;
s.t. c6: x[3] == 2;
