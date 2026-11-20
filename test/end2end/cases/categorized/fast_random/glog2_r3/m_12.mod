var x{1..4} integer >= 0 <= 3;

maximize o: 0*x[1] + 0*x[2] - 4*x[3] + x[4];

s.t. c1: atmost 3 {i in 1..4} (x[i] >= 0);
s.t. c2: alldiff(x[1], x[3], x[4], x[2]);
s.t. c3: (if x[4] == 2 then x[2] else x[1]) >= 0;
