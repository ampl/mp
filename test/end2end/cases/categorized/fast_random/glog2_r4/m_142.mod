var x{1..6} integer >= -2 <= 1;

minimize o: 5*x[1] + 0*x[2] + 0*x[3] - x[4] - 4*x[5] + 3*x[6];

s.t. c1: max(x[1], x[1]) - min(x[1], x[1]) < -2;
