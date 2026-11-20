var x{1..5} integer >= -3 <= 1;
minimize o: (-3)*x[1] + (-4)*x[2] + (-2)*x[3] + (0)*x[4] + (5)*x[5];
s.t. c1: atmost 2 {i in 1..5} (x[i] >= -2);
s.t. c2: max(x[2], x[2]) - min(x[2], x[1]) <= 1;
s.t. c3: (if x[1] == 1 then x[1] else x[4]) > -3;
s.t. c4: atleast 3 {i in 1..5} (x[i] >= 0);
