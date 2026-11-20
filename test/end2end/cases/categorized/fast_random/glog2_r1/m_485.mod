var x{1..5} integer >= 0 <= 4;
maximize o: (-3)*x[1] + (-4)*x[2] + (5)*x[3] + (-2)*x[4] + (4)*x[5];
s.t. c1: (x[4] != 4) ==> (x[4] == 1);
s.t. c2: max(x[4], x[5]) - min(x[4], x[2]) == 1;
s.t. c3: (x[3] <= 1) ==> (x[5] == 0) else (x[1] >= 2);
s.t. c4: numberof 4 in (x[1], x[2], x[3], x[4], x[5]) != 2;
s.t. c5: atmost 3 {i in 1..5} (x[i] >= 4);
