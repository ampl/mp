var x{1..4} integer >= -2 <= 2;

maximize o: -5*x[1] - 2*x[2] + 2*x[3] - 3*x[4];

s.t. c1: atleast 2 {i in 1..4} (x[i] >= 2);
