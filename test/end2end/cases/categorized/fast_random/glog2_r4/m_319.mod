var x{1..5} integer >= 0 <= 4;

maximize o: 3*x[1] - 2*x[2] - 4*x[3] + x[4] + 0*x[5];

s.t. c1: not ((x[1] != 2) and (x[2] != 4));
s.t. c2: (x[2] <= 3) ==> (x[2] >= 0) else (x[4] >= 2);
s.t. c3: atmost 3 {i in 1..5} (x[i] >= 3);
s.t. c4: not ((x[1] != 2) and (x[2] != 4));
s.t. c5: atmost 3 {i in 1..5} (x[i] >= 3);
