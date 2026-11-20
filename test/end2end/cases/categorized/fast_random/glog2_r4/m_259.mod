var x{1..4} integer >= 0 <= 4;
maximize o: (6)*x[1] + (4)*x[2] + (2)*x[3] + (-3)*x[4];
s.t. c1: atmost 3 {i in 1..4} (x[i] >= 4);
s.t. c2: numberof 1 in (x[1], x[2], x[3], x[4]) == 1;
s.t. c3: atleast 2 {i in 1..4} (x[i] >= 2);
s.t. c4: atmost 3 {i in 1..4} (x[i] >= 4);
s.t. c5: atleast 2 {i in 1..4} (x[i] >= 2);
