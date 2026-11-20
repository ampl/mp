var x{1..3} integer >= 0 <= 2;
minimize o: (-1)*x[1] + (2)*x[2] + (-4)*x[3];
s.t. c1: (x[3] < 2) ==> (x[3] != 2) else (x[3] > 2);
s.t. c2: (x[2] <= 2) or (x[1] == 2) or (x[2] >= 2);
s.t. c3: exists {i in 1..3} x[i] == 0;
s.t. c4: atmost 2 {i in 1..3} (x[i] >= 2);
s.t. c5: exists {i in 1..3} x[i] == 2;
