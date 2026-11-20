var x{1..6} integer >= -2 <= 2;

maximize o: 2*x[1] - 2*x[2] - 5*x[3] - 4*x[4] + 3*x[5] - 3*x[6];

s.t. c1: atleast 3 {i in 1..6} (x[i] >= -2);
