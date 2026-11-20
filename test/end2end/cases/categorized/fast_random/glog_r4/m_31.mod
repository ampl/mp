var x{1..3} integer >= 0 <= 2;
maximize o: (-2)*x[1] + (4)*x[2] + (-3)*x[3];
s.t. c1: atmost 1 {i in 1..3} (x[i] >= 1);
s.t. c2: max(x[1], x[1]) - min(x[2], x[2]) != 1;
s.t. c3: forall {i in 1..2} x[i] != x[i+1];
