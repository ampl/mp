var x{1..4} integer >= 0 <= 2;
maximize o: (6)*x[1] + (-3)*x[2] + (5)*x[3] + (-3)*x[4];
s.t. c1: atmost 1 {i in 1..4} (x[i] >= 1);
s.t. c2: max(x[1], x[2]) - min(x[4], x[1]) <= 0;
s.t. c3: abs(x[3] - x[4]) == 1;
s.t. c4: forall {i in 1..3} x[i] != x[i+1];
s.t. c5: atmost 2 {i in 1..4} (x[i] >= 2);
s.t. c6: ((x[2] != 2) ==> (x[3] != 1)) and ((x[3] != 1) ==> (x[4] >= 2));
s.t. c7: atleast 3 {i in 1..4} (x[i] >= 1);
