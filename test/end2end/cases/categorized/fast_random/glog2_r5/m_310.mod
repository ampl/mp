var x{1..4} integer >= 0 <= 3;

maximize o: -2*x[1] - 4*x[2] + 0*x[3] - 4*x[4];

s.t. c1: atleast 3 {i in 1..4} (x[i] >= 1);
