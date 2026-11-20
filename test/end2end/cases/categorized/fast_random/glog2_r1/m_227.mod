var x{1..4} integer >= 0 <= 3;

maximize o: 17*x[1] - 2000*x[2] + 3000*x[3] - 333*x[4];

s.t. c1: atmost 2 {i in 1..4} (x[i] >= 3);
