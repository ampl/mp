var x{1..5} integer >= 0 <= 4;

minimize o: -3*x[1] + 6*x[2] + 3*x[3] + 3*x[4] + 0*x[5];

s.t. c1: atleast 2 {i in 1..5} (x[i] >= 4);
