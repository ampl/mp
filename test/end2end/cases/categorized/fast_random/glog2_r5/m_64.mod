var x{1..4} integer >= -3 <= 2;

maximize o: -4*x[1] - 2*x[2] - x[3] + 0*x[4];

s.t. c1: atmost 3 {i in 1..4} (x[i] >= 1);
