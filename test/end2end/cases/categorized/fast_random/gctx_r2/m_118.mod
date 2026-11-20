var x{1..4} integer >= 0 <= 3;

minimize o: -3*x[1] + 0*x[2] + 2*x[3] - 3*x[4];

s.t. c1: not (x[1] <= x[1]);
s.t. c2: not (x[1] <= x[1]);
s.t. c3: not ((x[2] > x[3]) and (x[1] <= x[1]));
s.t. c4: not (max(x[1], x[3]) == 2);
