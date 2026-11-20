var x{1..3} integer >= -2 <= 2;

minimize o: -4*x[1] + x[2] - 3*x[3];

s.t. c1: abs(x[1] - x[2]) > 3;
s.t. c2: x[1] == -2;
s.t. c3: x[2] == 2;
