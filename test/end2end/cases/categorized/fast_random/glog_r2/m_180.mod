var x{1..3} integer >= -2 <= 1;

minimize o: 3*x[1] + 5*x[2] + 4*x[3];

s.t. c1: (x[2] > -2) or (x[1] <= 0) or (x[1] == 1);
s.t. c2: numberof -2 in (x[1], x[2], x[3]) < 3;
