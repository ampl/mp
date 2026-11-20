var x{1..3} integer >= 0 <= 2;

minimize o: 3*x[1] + 2*x[2] - 4*x[3];

s.t. c1: max(x[1], x[2]) - min(x[1], x[3]) == 1;
s.t. c2: abs(x[1] - x[1]) > 0;
