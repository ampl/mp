var x{1..4} integer >= -3 <= 1;
minimize o: (2)*x[1] + (5)*x[2] + (1)*x[3] + (-4)*x[4];
s.t. c1: (x[2] >= -1) <==> (x[3] + x[3] == -1);
s.t. c2: atmost 1 {i in 1..4} (x[i] >= 1);
s.t. c3: (x[1] < 0) and (x[1] == 1) ==> (x[4] >= 0);
s.t. c4: exists {i in 1..4} x[i] == 1;
