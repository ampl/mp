var x{1..3} integer >= 0 <= 4;

minimize o: -4*x[1] + 2*x[2] + 3*x[3];

s.t. c1: max(x[1], x[2]) - min(x[2], x[2]) <= 4;
