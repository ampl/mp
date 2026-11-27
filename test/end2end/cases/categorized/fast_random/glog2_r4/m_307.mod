var x{1..6} integer >= 0 <= 2;

maximize o: 1000*x[1] + 2000*x[2] - 666*x[3] + 3000*x[4] + 51*x[5] + 17*x[6];

s.t. c1: atmost 3 {i in 1..6} (x[i] >= 2);
