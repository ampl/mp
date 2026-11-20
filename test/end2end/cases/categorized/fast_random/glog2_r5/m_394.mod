var x{1..5} integer >= -3 <= 2;

minimize o: -x[1] - x[2] + 2*x[3] - 3*x[4] - 5*x[5];

s.t. c1: atleast 2 {i in 1..5} (x[i] >= 0);
