var x{1..4} integer >= 0 <= 3;

maximize o: -x[1] + 3*x[2] - 4*x[3] - 2*x[4];

s.t. c1: atleast 3 {i in 1..4} (x[i] >= 0);
