var x{1..5} integer >= -3 <= 1;

minimize o: 3*x[1] + 4*x[2] + 5*x[3] - 3*x[4] + x[5];

s.t. c1: abs(x[4] - x[2]) >= 0;
s.t. c2: numberof -1 in (x[1], x[2], x[3], x[4], x[5]) < 3;
s.t. c3: max(x[2], x[2]) - min(x[5], x[1]) < -1;
