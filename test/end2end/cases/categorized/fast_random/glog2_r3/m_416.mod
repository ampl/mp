var x{1..5} integer >= -2 <= 0;

maximize o: 0*x[1] + 0*x[2] + 0*x[3] + 0*x[4] + 0*x[5];

s.t. c1: atmost 2 {i in 1..5} (x[i] >= -1);
