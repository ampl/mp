var x{1..4} integer >= 0 <= 4;
minimize o: (-3)*x[1] + (-3)*x[2] + (-2)*x[3] + (-1)*x[4];
s.t. c1: numberof 2 in (x[1], x[2], x[3], x[4]) < 3;
s.t. c2: ((x[2] < 4) ==> (x[1] >= 3)) and ((x[1] >= 3) ==> (x[1] >= 4));
s.t. c3: not ((x[3] > 0) and (x[4] > 0));
s.t. c4: (x[2] == 2) ==> (x[1] <= 2) else (x[3] > 2);
s.t. c5: forall {i in 1..3} x[i] != x[i+1];
