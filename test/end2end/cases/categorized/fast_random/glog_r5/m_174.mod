var x{1..3} integer >= 0 <= 2;
minimize o: (-5)*x[1] + (2)*x[2] + (-1)*x[3];
s.t. c1: ((if x[2] > 0 then 1 else 0) + (if x[1] >= 0 then 1 else 0)) >= 1;
s.t. c2: forall {i in 1..2} x[i] != x[i+1];
s.t. c3: (x[3] < 1) ==> (x[1] != 2);
s.t. c4: exists {i in 1..3} x[i] == 2;
