var x{1..4} integer >= -2 <= 0;

minimize o: -x[1] - 3*x[2] + 5*x[3] + 5*x[4];

s.t. c1: abs(x[1] - x[4]) <= 2;
s.t. c2: (if x[1] != 0 then x[1] else x[2]) >= -1;
s.t. c3: (x[3] == 0) or (x[3] < 0) or (x[3] != 0);
