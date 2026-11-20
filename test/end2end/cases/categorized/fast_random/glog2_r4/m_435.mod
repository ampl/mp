var x{1..5} integer >= -2 <= 2;
maximize o: (-4)*x[1] + (3)*x[2] + (1)*x[3] + (-3)*x[4] + (0)*x[5];
s.t. c1: atmost 3 {i in 1..5} (x[i] >= -1);
s.t. c2: (x[2] <= -1) or (x[2] <= 1) or (x[4] > -2);
s.t. c3: abs(x[4] - x[4]) != 3;
s.t. c4: max(x[2], x[1]) - min(x[5], x[1]) != -2;
s.t. c5: (if x[2] < -1 then x[5] else x[1]) < 0;
s.t. c6: numberof -1 in (x[1], x[2], x[3], x[4], x[5]) <= 3;
