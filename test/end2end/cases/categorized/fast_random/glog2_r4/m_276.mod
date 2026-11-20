var x{1..5} integer >= 0 <= 4;

maximize o: 0*x[1] + 0*x[2] + 0*x[3] + 0*x[4] + 0*x[5];

s.t. c1: (x[3] >= 2) ==> (x[1] == 1);
s.t. c2: (x[5] != 4) <==> (x[3] + x[4] > 4);
s.t. c3: abs(x[4] - x[1]) == 3;
s.t. c4: numberof 1 in (x[1], x[2], x[3], x[4], x[5]) == 2;
