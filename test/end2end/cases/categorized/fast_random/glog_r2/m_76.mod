var x{1..5} integer >= 0 <= 2;
minimize o: (-3)*x[1] + (-5)*x[2] + (3)*x[3] + (2)*x[4] + (-4)*x[5];
s.t. c1: not ((x[2] >= 0) and (x[1] != 0));
s.t. c2: forall {i in 1..4} x[i] != x[i+1];
s.t. c3: (x[3] != 1) and (x[4] >= 2) ==> (x[2] >= 1);
