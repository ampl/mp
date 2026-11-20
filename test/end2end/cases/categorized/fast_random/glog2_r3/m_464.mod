var x{1..4} integer >= -3 <= 1;

maximize o: -2*x[1] - 2*x[2] - x[3] + 5*x[4];

s.t. c1: atleast 2 {i in 1..4} (x[i] >= 0);
