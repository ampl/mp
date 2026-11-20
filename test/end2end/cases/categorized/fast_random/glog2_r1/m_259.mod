var x{1..5} integer >= -2 <= 0;
minimize o: (0)*x[1] + (4)*x[2] + (-3)*x[3] + (-4)*x[4] + (-5)*x[5];
s.t. c1: forall {i in 1..4} x[i] != x[i+1];
s.t. c2: exists {i in 1..5} x[i] == -1;
s.t. c3: max(x[3], x[3]) - min(x[2], x[3]) != -2;
s.t. c4: forall {i in 1..4} x[i] != x[i+1];
s.t. c5: max(x[3], x[3]) - min(x[2], x[3]) != -2;
