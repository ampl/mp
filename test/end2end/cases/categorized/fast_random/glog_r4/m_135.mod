var x{1..4} integer >= -2 <= 0;
minimize o: (4)*x[1] + (4)*x[2] + (-2)*x[3] + (2)*x[4];
s.t. c1: max(x[1], x[4]) - min(x[3], x[3]) >= -2;
s.t. c2: max(x[2], x[4]) - min(x[3], x[1]) < -2;
s.t. c3: max(x[2], x[3]) - min(x[4], x[3]) == -1;
s.t. c4: (if x[1] >= -1 then x[2] else x[4]) <= -2;
