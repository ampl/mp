var x{1..5} integer >= 0 <= 3;

minimize o: 3*x[1] + x[2] - 5*x[3] + 5*x[4] - 5*x[5];

s.t. c1: (x[4] >= 0) or (x[4] == 3) or (x[2] <= 3);
s.t. c2: numberof 1 in (x[1], x[2], x[3], x[4], x[5]) == 1;
