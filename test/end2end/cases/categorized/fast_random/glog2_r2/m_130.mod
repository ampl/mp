var x{1..5} integer >= 0 <= 4;
maximize o: (-4)*x[1] + (-5)*x[2] + (3)*x[3] + (-4)*x[4] + (1)*x[5];
s.t. c1: forall {i in 1..4} x[i] != x[i+1];
s.t. c2: (if x[4] < 4 then x[4] else x[1]) < 1;
s.t. c3: (x[2] > 0) or (x[5] > 0) or (x[2] != 4);
s.t. c4: exists {i in 1..5} x[i] == 3;
