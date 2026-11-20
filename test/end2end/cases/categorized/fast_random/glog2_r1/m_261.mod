var x{1..5} integer >= -2 <= 1;
minimize o: (-4)*x[1] + (2)*x[2] + (3)*x[3] + (-3)*x[4] + (-3)*x[5];
s.t. c1: forall {i in 1..4} x[i] != x[i+1];
s.t. c2: max(x[2], x[5]) - min(x[3], x[5]) == -3;
s.t. c3: max(x[5], x[3]) - min(x[2], x[5]) > -3;
