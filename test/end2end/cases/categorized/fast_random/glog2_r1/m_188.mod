var x{1..5} integer >= 0 <= 1;
maximize o: (3)*x[1] + (-3)*x[2] + (3)*x[3] + (-2)*x[4] + (0)*x[5];
s.t. c1: atmost 1 {i in 1..5} (x[i] >= 1);
s.t. c2: ((if x[1] >= 1 then 1 else 0) + (if x[2] >= 1 then 1 else 0)) >= 1;
