var x{1..5} integer >= 0 <= 4;

maximize o: 4*x[1] + 0*x[2] - 3*x[3] + 2*x[4] + 3*x[5];

s.t. c1: (if x[2] != 0 then x[1] else x[4]) >= 4;
s.t. c2: atmost 2 {i in 1..5} (x[i] >= 0);
