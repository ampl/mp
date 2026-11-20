var x{1..5} integer >= -2 <= 1;

maximize o: x[1] - 5*x[2] - 5*x[3] - 5*x[4] + x[5];

s.t. c1: alldiff(x[1], x[2], x[4]);
s.t. c2: (x[5] >= 2) ==> (x[3] != -3);
s.t. c3: atmost 1 {i in 1..5} (x[i] >= -2);
