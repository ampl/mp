var x{1..6} integer >= 0 <= 4;

minimize o: 0*x[1] - 4*x[2] - 3*x[3] - x[4] + 2*x[5] - 4*x[6];

s.t. c1: abs(x[3] - x[3]) > 2;
