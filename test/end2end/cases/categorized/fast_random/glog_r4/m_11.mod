var x{1..4} integer >= -2 <= 1;

minimize o: 0*x[1] - x[2] + x[3] + 6*x[4];

s.t. c1: atleast 3 {i in 1..4} (x[i] >= 0);
