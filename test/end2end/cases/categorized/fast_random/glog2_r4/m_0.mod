var x{1..5} integer >= 0 <= 3;
minimize o: (-3)*x[1] + (-5)*x[2] + (-1)*x[3] + (2)*x[4] + (5)*x[5];
s.t. c1: forall {i in 1..4} x[i] != x[i+1];
s.t. c2: alldiff(x[1], x[2], x[5], x[3]);
