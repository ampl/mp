var x{1..4} integer >= -2 <= 1;

minimize o: 6*x[1] + x[2] - 3*x[3] - 5*x[4];

s.t. c1: numberof -2 in (x[1], x[2], x[3], x[4]) <= 2;
s.t. c2: (x[3] > -2) and (x[3] <= 0) ==> (x[4] >= -2);
s.t. c3: abs(x[4] - x[2]) < 2;
