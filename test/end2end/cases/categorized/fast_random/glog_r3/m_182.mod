var x{1..4} integer >= 0 <= 4;

maximize o: 2*x[1] + 6*x[2] - 4*x[3] + 2*x[4];

s.t. c1: max(x[2], x[2]) - min(x[2], x[1]) <= 4;
s.t. c2: (x[2] >= 3) or (x[4] != 2) or (x[3] >= 1);
s.t. c3: atmost 2 {i in 1..4} (x[i] >= 4);
s.t. c4: ((if x[2] > 2 then 1 else 0) + (if x[2] >= 2 then 1 else 0)) >= 1;
