var x{1..4} integer >= -2 <= 2;
minimize o: (-1)*x[1] + (4)*x[2] + (1)*x[3] + (1)*x[4];
s.t. c1: numberof -1 in (x[1], x[2], x[3], x[4]) >= 3;
s.t. c2: exists {i in 1..4} x[i] == 0;
s.t. c3: atleast 3 {i in 1..4} (x[i] >= 2);
s.t. c4: atmost 2 {i in 1..4} (x[i] >= -2);
