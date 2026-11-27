var x{1..5} integer >= -2 <= 2;

maximize o: 0*x[1] + 3*x[2] + 3*x[3] + x[4] + 5*x[5];

s.t. c1: atleast 1 {i in 1..5} (x[i] >= -1);
