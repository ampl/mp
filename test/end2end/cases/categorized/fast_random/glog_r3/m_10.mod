var x{1..5} integer >= -2 <= 1;
minimize o: (-3)*x[1] + (0)*x[2] + (-5)*x[3] + (1)*x[4] + (0)*x[5];
s.t. c1: numberof 1 in (x[1], x[2], x[3], x[4], x[5]) == 1;
s.t. c2: abs(x[5] - x[4]) < 1;
s.t. c3: forall {i in 1..4} x[i] != x[i+1];
s.t. c4: (x[1] >= -2) ==> (x[1] < 1);
