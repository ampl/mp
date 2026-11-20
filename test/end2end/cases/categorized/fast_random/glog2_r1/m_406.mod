var x{1..3} integer >= -2 <= 2;

maximize o: 0*x[1] + 0*x[2] + 0*x[3];

s.t. c1: numberof 2 in (x[1], x[2], x[3]) == 0;
s.t. c2: (x[3] > -2) <==> (x[1] + x[3] != 1);
