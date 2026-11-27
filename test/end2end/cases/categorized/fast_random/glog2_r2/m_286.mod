var x{1..4} integer >= 0 <= 2;

maximize o: 0*x[1] + 0*x[2] + 0*x[3] + 0*x[4];

s.t. c1: (x[4] != 0) or (x[3] <= 0) or (x[4] == 0);
s.t. c2: (x[4] <= 2) <==> (x[3] + x[1] <= 2);
s.t. c3: atmost 2 {i in 1..4} (x[i] >= 2);
s.t. c4: not ((x[2] <= 2) and (x[3] != 1));
