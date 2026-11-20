var x{1..3} integer >= -2 <= 2;
maximize o: (-4)*x[1] + (-3)*x[2] + (-4)*x[3];
s.t. c1: ((x[3] == 0) ==> (x[1] == 2)) and ((x[1] == 2) ==> (x[3] >= 0));
s.t. c2: max(x[3], x[2]) - min(x[1], x[1]) == -2;
s.t. c3: (if x[3] == 1 then x[3] else x[2]) > -1;
s.t. c4: atmost 1 {i in 1..3} (x[i] >= -1);
