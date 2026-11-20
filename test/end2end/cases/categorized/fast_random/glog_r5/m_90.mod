var x{1..6} integer >= 0 <= 2;

minimize o: 4*x[1] + 5*x[2] + x[3] + 6*x[4] + 0*x[5] - 2*x[6];

s.t. c1: atleast 2 {i in 1..6} (x[i] >= 0);
