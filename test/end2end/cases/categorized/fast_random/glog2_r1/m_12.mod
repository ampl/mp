var x{1..3} integer >= 0 <= 4;
minimize o: (2)*x[1] + (6)*x[2] + (-5)*x[3];
s.t. c1: exists {i in 1..3} x[i] == 0;
s.t. c2: max(x[1], x[1]) - min(x[1], x[1]) >= 3;
s.t. c3: atmost 2 {i in 1..3} (x[i] >= 2);
