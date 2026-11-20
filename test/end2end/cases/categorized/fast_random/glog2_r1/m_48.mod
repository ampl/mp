var x{1..7} integer >= 0 <= 1;
maximize o: (3)*x[1] + (0)*x[2] + (-1)*x[3] + (2)*x[4] + (-2)*x[5] + (6)*x[6] + (1)*x[7];
s.t. c1: abs(x[4] - x[7]) == 0;
s.t. c2: (if x[7] != 1 then x[7] else x[4]) <= 0;
s.t. c3: atleast 3 {i in 1..7} (x[i] >= 0);
s.t. c4: atmost 2 {i in 1..7} (x[i] >= 1);
