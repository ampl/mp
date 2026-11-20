var x{1..5} integer >= 0 <= 1;
minimize o: (2)*x[1] + (3)*x[2] + (2)*x[3] + (0)*x[4] + (0)*x[5];
s.t. c1: atmost 2 {i in 1..5} (x[i] >= 0);
s.t. c2: max(x[4], x[4]) - min(x[3], x[4]) <= 1;
s.t. c3: (x[4] > 1) <==> (x[2] + x[4] <= 0);
