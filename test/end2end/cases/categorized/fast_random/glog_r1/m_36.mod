var x{1..4} integer >= -2 <= 2;

maximize o: -3*x[1] - x[2] + 0*x[3] + 2*x[4];

s.t. c1: (x[2] == -2) or (x[4] >= 2) or (x[1] <= 2);
s.t. c2: atmost 3 {i in 1..4} (x[i] >= 2);
