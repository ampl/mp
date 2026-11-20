var x{1..5} integer >= 0 <= 2;
minimize o: (0)*x[1] + (-3)*x[2] + (-3)*x[3] + (3)*x[4] + (1)*x[5];
s.t. c1: (x[1] == 0) <==> (x[1] + x[4] > 2);
s.t. c2: max(x[2], x[2]) - min(x[3], x[2]) < 2;
s.t. c3: atmost 1 {i in 1..5} (x[i] >= 0);
s.t. c4: atmost 1 {i in 1..5} (x[i] >= 2);
