var x{1..3} integer >= -3 <= 1;
minimize o: (2)*x[1] + (-5)*x[2] + (3)*x[3];
s.t. c1: forall {i in 1..2} x[i] != x[i+1];
s.t. c2: (x[2] == -3) ==> (x[3] != -1) else (x[3] >= -1);
s.t. c3: forall {i in 1..2} x[i] != x[i+1];
s.t. c4: numberof -1 in (x[1], x[2], x[3]) < 0;
