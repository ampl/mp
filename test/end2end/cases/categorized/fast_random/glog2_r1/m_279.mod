var x{1..4} integer >= 0 <= 2;
minimize o: (-5)*x[1] + (0)*x[2] + (-5)*x[3] + (1)*x[4];
s.t. c1: exists {i in 1..4} x[i] == 0;
s.t. c2: (x[2] <= 2) and (x[2] == 0) ==> (x[2] >= 2);
s.t. c3: forall {i in 1..3} x[i] != x[i+1];
s.t. c4: (x[2] != 0) or (x[2] <= 2) or (x[1] == 0);
s.t. c5: exists {i in 1..4} x[i] == 0;
s.t. c6: (x[2] != 0) or (x[2] <= 2) or (x[1] == 0);
