var x{1..5} integer >= -3 <= 2;

minimize o: -3*x[1] + 0*x[2] - 3*x[3] - 2*x[4] - 4*x[5];

s.t. c1: numberof -1 in (x[1], x[2], x[3], x[4], x[5]) >= 1;
s.t. c2: (if x[4] < -1 then x[3] else x[3]) < -1;
s.t. c3: not ((x[1] >= 0) and (x[3] == 2));
