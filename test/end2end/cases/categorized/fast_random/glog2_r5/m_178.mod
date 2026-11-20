var x{1..5} integer >= 0 <= 1;

maximize o: x[1] + 3*x[2] + 4*x[3] + x[4] - 2*x[5];

s.t. c1: max(x[4], x[2]) - min(x[5], x[1]) == 0;
s.t. c2: numberof 0 in (x[1], x[2], x[3], x[4], x[5]) > 2;
s.t. c3: (x[2] != 0) or (x[2] > 0) or (x[2] != 1);
