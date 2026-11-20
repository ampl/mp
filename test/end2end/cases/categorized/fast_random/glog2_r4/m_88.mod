var x{1..7} integer >= 0 <= 1;

maximize o: 2*x[1] + 6*x[2] + 2*x[3] - 5*x[4] + 2*x[5] + 3*x[6] + 3*x[7];

s.t. c1: exists {i in 1..7} x[i] == 1;
s.t. c2: atmost 2 {i in 1..7} (x[i] >= 1);
