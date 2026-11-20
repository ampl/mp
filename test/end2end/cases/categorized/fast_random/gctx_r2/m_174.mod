var x{1..4} integer >= -1 <= 2;

minimize o: x[1] + 6*x[2] + 2*x[3] - 3*x[4];

s.t. c1: not ((abs(x[4] - x[4]) <= 0) and (abs(x[4] - x[4]) <= 0));
s.t. c2: not (abs(x[4] - x[4]) <= 0);
s.t. c3: not ((x[3] >= 1) and (x[3] >= 1));
