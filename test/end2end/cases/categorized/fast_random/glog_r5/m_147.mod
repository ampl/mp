var x{1..3} integer >= -2 <= 1;
minimize o: (-5)*x[1] + (2)*x[2] + (-2)*x[3];
s.t. c1: exists {i in 1..3} x[i] == -2;
s.t. c2: forall {i in 1..2} x[i] != x[i+1];
s.t. c3: (x[3] == -2) and (x[2] < 0) ==> (x[3] >= -2);
