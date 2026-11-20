var x{1..7} integer >= 0 <= 1;

maximize o: -4*x[1] - x[2] + 6*x[3] - x[4] - x[5] - 5*x[6] - 4*x[7];

s.t. c1: numberof 1 in (x[1], x[2], x[3], x[4], x[5], x[6], x[7]) >= 0;
s.t. c2: (x[1] != 1) <==> (x[5] + x[5] == 0);
