var x{1..5} integer >= 0 <= 2;

maximize o: -2*x[1] - 2*x[2] + 6*x[3] - 4*x[4] + 2*x[5];

s.t. c1: atmost 2 {i in 1..5} (x[i] >= 2);
s.t. c2: not ((x[3] >= 1) and (x[1] == 1));
