var x{1..7} integer >= 0 <= 1;

maximize o: 6*x[1] + 2*x[2] + 4*x[3] + x[4] - 4*x[5] - 5*x[6] - 5*x[7];

s.t. c1: exists {i in 1..7} x[i] == 0;
s.t. c2: numberof 0 in (x[1], x[2], x[3], x[4], x[5], x[6], x[7]) != 0;
