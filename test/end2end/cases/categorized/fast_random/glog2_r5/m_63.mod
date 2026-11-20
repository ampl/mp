var x{1..5} integer >= -2 <= 1;
maximize o: (1)*x[1] + (-3)*x[2] + (-2)*x[3] + (-1)*x[4] + (-5)*x[5];
s.t. c1: not ((x[4] > 1) and (x[3] < -2));
s.t. c2: forall {i in 1..4} x[i] != x[i+1];
s.t. c3: exists {i in 1..5} x[i] == 1;
s.t. c4: x[4] == -1;
s.t. c5: x[5] == -1;
