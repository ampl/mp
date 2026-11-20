var x{1..3} integer >= 0 <= 4;
minimize o: (-4)*x[1] + (-1)*x[2] + (1)*x[3];
s.t. c1: exists {i in 1..3} x[i] == 2;
s.t. c2: max(x[1], x[1]) - min(x[2], x[1]) == 4;
s.t. c3: (x[3] == 2) ==> (x[2] >= 2) else (x[1] != 2);
s.t. c4: (x[1] == 2) or (x[3] > 2) or (x[1] <= 0);
