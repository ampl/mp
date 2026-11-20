var x{1..6} integer >= -2 <= 1;

minimize o: -4*x[1] - x[2] + 2*x[3] + 5*x[4] + 4*x[5] + 0*x[6];

s.t. c1: atleast 2 {i in 1..6} (x[i] >= -2);
s.t. c2: atleast 2 {i in 1..6} (x[i] >= -2);
