var x{1..5} integer >= -2 <= 3;
maximize o: (-2)*x[1] + (4)*x[2] + (0)*x[3] + (3)*x[4] + (-1)*x[5];
s.t. c1: numberof -2 in (x[1], x[2], x[3], x[4], x[5]) < 0;
s.t. c2: atmost 1 {i in 1..5} (x[i] >= 0);
s.t. c3: abs(x[2] - x[4]) > 2;
s.t. c4: (if x[3] >= -2 then x[5] else x[3]) >= 0;
