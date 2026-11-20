var x{1..3} integer >= 0 <= 4;
maximize o: (-3)*x[1] + (-4)*x[2] + (4)*x[3];
s.t. c1: atmost 3 {i in 1..3} (x[i] >= 4);
