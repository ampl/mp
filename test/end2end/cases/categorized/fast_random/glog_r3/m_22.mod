var x{1..6} integer >= 0 <= 3;

minimize o: x[1] + 6*x[2] - 4*x[3] + 6*x[4] + 5*x[5] - 4*x[6];

s.t. c1: max(x[2], x[1]) - min(x[2], x[4]) < 0;
s.t. c2: (x[6] > 2) or (x[5] < 1) or (x[3] <= 3);
s.t. c3: numberof 2 in (x[1], x[2], x[3], x[4], x[5], x[6]) == 2;
