var x{1..4} integer >= 0 <= 2;

maximize o: 0*x[1] + 0*x[2] + 0*x[3] + 0*x[4];

s.t. c1: numberof 0 in (x[1], x[2], x[3], x[4]) < 2;
s.t. c2: (x[4] <= 1) <==> (x[3] + x[1] < 0);
