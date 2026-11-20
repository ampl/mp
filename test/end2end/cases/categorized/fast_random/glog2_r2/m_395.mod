var x{1..4} integer >= -2 <= 0;
maximize o: (2)*x[1] + (3)*x[2] + (1)*x[3] + (2)*x[4];
s.t. c1: numberof -1 in (x[1], x[2], x[3], x[4]) != 3;
s.t. c2: (x[3] != 0) or (x[2] == 0) or (x[4] <= 0);
s.t. c3: max(x[3], x[4]) - min(x[4], x[3]) >= -2;
s.t. c4: (if x[3] < -2 then x[4] else x[4]) >= 0;
s.t. c5: ((if x[4] == -1 then 1 else 0) + (if x[2] >= -1 then 1 else 0)) >= 1;
s.t. c6: forall {i in 1..3} x[i] != x[i+1];
s.t. c7: atmost 2 {i in 1..4} (x[i] >= -2);
