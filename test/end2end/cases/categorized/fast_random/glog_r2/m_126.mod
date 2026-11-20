var x{1..6} integer >= -2 <= 0;

minimize o: -4*x[1] - 2*x[2] + 0*x[3] + 6*x[4] - 2*x[5] + 5*x[6];

s.t. c1: exists {i in 1..6} x[i] == 0;
s.t. c2: abs(x[6] - x[3]) < 2;
s.t. c3: abs(x[1] - x[5]) > 2;
s.t. c4: (x[4] >= -2) or (x[5] > 0) or (x[4] != -2);
