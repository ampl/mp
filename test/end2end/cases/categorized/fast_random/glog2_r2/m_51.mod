var x{1..4} integer >= 0 <= 3;
minimize o: (-2)*x[1] + (-4)*x[2] + (5)*x[3] + (-5)*x[4];
s.t. c1: forall {i in 1..3} x[i] != x[i+1];
s.t. c2: (x[3] > -1) ==> (x[2] >= 3);
s.t. c3: atmost 2 {i in 1..4} (x[i] >= -1);
s.t. c4: not ((x[2] >= -1) and (x[3] != 0));
