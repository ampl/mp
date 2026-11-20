var x{1..5} integer >= 0 <= 1;
minimize o: (-2)*x[1] + (1)*x[2] + (-3)*x[3] + (3)*x[4] + (-5)*x[5];
s.t. c1: (x[5] >= 1) ==> (x[1] >= 0) else (x[1] >= 1);
s.t. c2: numberof 0 in (x[1], x[2], x[3], x[4], x[5]) <= 1;
s.t. c3: max(x[1], x[1]) - min(x[1], x[1]) == 0;
s.t. c4: atleast 3 {i in 1..5} (x[i] >= 1);
