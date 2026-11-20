var x{1..4} integer >= -3 <= 1;

minimize o: -3*x[1] + 0*x[2] - 5*x[3] - 2*x[4];

s.t. c1: atleast 3 {i in 1..4} (x[i] >= -3);
