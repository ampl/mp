var x{1..6} integer >= 0 <= 2;

minimize o: 2*x[1] + 2*x[2] + 2*x[3] + 5*x[4] + 3*x[5] - x[6];

s.t. c1: exists {i in 1..6} x[i] == 1;
s.t. c2: abs(x[3] - x[3]) < 0;
s.t. c3: abs(x[3] - x[2]) <= 0;
