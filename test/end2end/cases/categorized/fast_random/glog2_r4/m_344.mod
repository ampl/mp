var x{1..5} integer >= -2 <= 3;

minimize o: 3*x[1] + 6*x[2] + 6*x[3] + 6*x[4] - x[5];

s.t. c1: (x[2] != -2) or (x[5] <= -2) or (x[3] >= 2);
