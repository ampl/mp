var x{1..5} integer >= -2 <= 2;

maximize o: 2*x[1] - 3*x[2] + x[3] + 6*x[4] + 3*x[5];

s.t. c1: numberof -1 in (x[1], x[2], x[3], x[4], x[5]) == 3;
s.t. c2: (x[3] == 1) or (x[5] == -1) or (x[4] >= -2);
