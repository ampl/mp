var x{1..5} integer >= -2 <= 3;
minimize o: (-5)*x[1] + (4)*x[2] + (2)*x[3] + (2)*x[4] + (3)*x[5];
s.t. c1: (if x[4] <= -2 then x[1] else x[5]) == 2;
s.t. c2: (x[5] > 0) <==> (x[4] + x[5] >= 1);
s.t. c3: atmost 3 {i in 1..5} (x[i] >= 3);
