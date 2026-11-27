var x{1..3} integer >= 0 <= 2;

minimize o: 4*x[1] - 2*x[2] + 5*x[3];

s.t. c1: (x[3] > 0) or (x[1] == 2) or (x[3] >= 1);
s.t. c2: not ((x[2] == 2) and (x[2] != 2));
s.t. c3: numberof 0 in (x[1], x[2], x[3]) <= 2;
