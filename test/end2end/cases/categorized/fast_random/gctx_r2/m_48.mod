var x{1..4} integer >= 0 <= 3;

minimize o: -x[1] + 3*x[2] + 0*x[3] + 0*x[4];

s.t. c1: (max(x[4], x[2]) < 0) or (max(x[1], x[2]) <= 3);
s.t. c2: not (max(x[4], x[2]) < 0);
s.t. c3: (max(x[4], x[2]) < 0) or (max(x[4], x[2]) < 0);
