var x{1..3} integer >= -2 <= 2;

maximize o: -4*x[1] + x[2] + 4*x[3];

s.t. c1: atmost 2 {i in 1..3} (x[i] >= -1);
