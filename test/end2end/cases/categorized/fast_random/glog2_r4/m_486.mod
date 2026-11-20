var x{1..4} integer >= 0 <= 4;

maximize o: 0*x[1] + 0*x[2] + 0*x[3] + 0*x[4];

s.t. c1: (x[3] < 1) or (x[1] < 4) or (x[3] <= 0);
s.t. c2: atmost 1 {i in 1..4} (x[i] >= 4);
s.t. c3: not ((x[1] < 4) and (x[3] >= 4));
