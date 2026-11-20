var x{1..5} integer >= -2 <= 1;

maximize o: 5*x[1] - 3*x[2] + 3*x[3] - x[4] + 6*x[5];

s.t. c1: (x[5] == -3) and (x[4] <= 2) ==> (x[3] >= -3);
s.t. c2: atmost 2 {i in 1..5} (x[i] >= 1);
