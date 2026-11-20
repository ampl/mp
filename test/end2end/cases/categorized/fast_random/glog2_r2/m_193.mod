var x{1..3} integer >= 0 <= 2;

minimize o: 2*x[1] + 3*x[2] + 6*x[3];

s.t. c1: abs(x[1] - x[1]) != 0;
s.t. c2: x[1] == 2;
