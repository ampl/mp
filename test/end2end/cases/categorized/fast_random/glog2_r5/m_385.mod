var x{1..3} integer >= -2 <= 2;
maximize o: (6)*x[1] + (2)*x[2] + (4)*x[3];
s.t. c1: atmost 1 {i in 1..3} (x[i] >= -1);
s.t. c2: (if x[1] < -1 then x[3] else x[1]) > 2;
s.t. c3: abs(x[1] - x[2]) <= 3;
s.t. c4: (x[1] <= 1) and (x[2] <= 0) ==> (x[3] >= 1);
s.t. c5: (x[3] > 2) and (x[2] >= 0) ==> (x[2] >= 2);
s.t. c6: alldiff(x[2], x[3]);
s.t. c7: numberof 2 in (x[1], x[2], x[3]) <= 2;
