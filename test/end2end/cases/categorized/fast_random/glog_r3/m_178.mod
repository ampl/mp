var x{1..4} integer >= 0 <= 3;

maximize o: -2*x[1] + x[2] - 4*x[3] - 2*x[4];

s.t. c1: (x[2] < 2) <==> (x[2] + x[2] < 1);
s.t. c2: not ((x[2] != 3) and (x[3] > 2));
s.t. c3: atmost 3 {i in 1..4} (x[i] >= 2);
