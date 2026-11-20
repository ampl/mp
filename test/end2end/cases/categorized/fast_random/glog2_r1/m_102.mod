var x{1..3} integer >= 0 <= 3;
minimize o: (3)*x[1] + (0)*x[2] + (-2)*x[3];
s.t. c1: atmost 3 {i in 1..3} (x[i] >= 3);
s.t. c2: forall {i in 1..2} x[i] != x[i+1];
s.t. c3: ((x[2] < 3) ==> (x[2] >= 3)) and ((x[2] >= 3) ==> (x[1] >= 3));
