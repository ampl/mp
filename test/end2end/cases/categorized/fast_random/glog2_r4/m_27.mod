var x{1..4} integer >= 0 <= 3;
maximize o: (17)*x[1] + (-999)*x[2] + (-2000)*x[3] + (-1000)*x[4];
s.t. c1: atmost 3 {i in 1..4} (x[i] >= 3);
s.t. c2: (if x[1] >= 3 then x[3] else x[4]) < 1;
s.t. c3: (x[3] < 0) and (x[1] != 0) ==> (x[3] >= 0);
