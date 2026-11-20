var x{1..5} integer >= -3 <= 1;
minimize o: (2)*x[1] + (6)*x[2] + (2)*x[3] + (-4)*x[4] + (3)*x[5];
s.t. c1: forall {i in 1..4} x[i] != x[i+1];
s.t. c2: numberof 0 in (x[1], x[2], x[3], x[4], x[5]) < 3;
s.t. c3: (x[1] >= -2) and (x[2] != 0) ==> (x[3] >= -2);
