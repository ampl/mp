var x{1..3} integer >= 0 <= 3;
minimize o: (1)*x[1] + (5)*x[2] + (-5)*x[3];
s.t. c1: atmost 1 {i in 1..3} (x[i] >= 3);
s.t. c2: (if x[2] < 1 then x[3] else x[1]) >= 0;
s.t. c3: numberof 1 in (x[1], x[2], x[3]) < 1;
s.t. c4: (x[3] >= 0) ==> (x[1] != 3) else (x[3] < 0);
