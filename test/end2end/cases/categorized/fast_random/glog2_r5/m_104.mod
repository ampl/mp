var x{1..5} integer >= -3 <= 2;

maximize o: -2*x[1] + x[2] + 0*x[3] + 3*x[4] - 4*x[5];

s.t. c1: max(x[2], x[5]) - min(x[1], x[4]) < -3;
