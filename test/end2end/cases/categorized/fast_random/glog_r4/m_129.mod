var x{1..5} integer >= -2 <= 0;

maximize o: -4*x[1] + 4*x[2] - 4*x[3] + 4*x[4] - x[5];

s.t. c1: atmost 2 {i in 1..5} (x[i] >= -2);
