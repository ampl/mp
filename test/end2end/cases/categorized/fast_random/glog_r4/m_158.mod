var x{1..5} integer >= 0 <= 4;
maximize o: (5)*x[1] + (1)*x[2] + (-4)*x[3] + (1)*x[4] + (5)*x[5];
s.t. c1: numberof 4 in (x[1], x[2], x[3], x[4], x[5]) <= 0;
s.t. c2: atmost 1 {i in 1..5} (x[i] >= 4);
s.t. c3: (if x[1] != 1 then x[4] else x[3]) < 3;
s.t. c4: (x[2] <= 4) or (x[2] <= 4) or (x[3] == 0);
