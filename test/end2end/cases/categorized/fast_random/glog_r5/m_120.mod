var x{1..5} integer >= -2 <= 0;

minimize o: x[1] - 4*x[2] - 5*x[3] - 4*x[4] + 3*x[5];

s.t. c1: atleast 1 {i in 1..5} (x[i] >= -1);
