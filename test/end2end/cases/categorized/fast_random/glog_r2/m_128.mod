var x{1..5} integer >= 0 <= 3;

maximize o: 6*x[1] - 5*x[2] + 0*x[3] + x[4] + 0*x[5];

s.t. c1: (x[5] != 2) or (x[1] >= 3) or (x[5] < 1);
s.t. c2: numberof 2 in (x[1], x[2], x[3], x[4], x[5]) != 1;
