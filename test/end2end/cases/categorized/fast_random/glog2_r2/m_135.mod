var x{1..3} integer >= 0 <= 3;
maximize o: (-2)*x[1] + (3)*x[2] + (-3)*x[3];
s.t. c1: atmost 2 {i in 1..3} (x[i] >= 1);
s.t. c2: ((if x[1] > 0 then 1 else 0) + (if x[2] >= 0 then 1 else 0)) >= 1;
s.t. c3: alldiff(x[1], x[2], x[3]);
s.t. c4: (x[3] < 0) or (x[1] > 1) or (x[1] != 2);
s.t. c5: atleast 3 {i in 1..3} (x[i] >= 3);
s.t. c6: atmost 1 {i in 1..3} (x[i] >= 3);
