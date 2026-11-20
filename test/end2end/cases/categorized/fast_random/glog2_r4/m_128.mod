var x{1..7} integer >= 0 <= 1;
minimize o: (0)*x[1] + (-5)*x[2] + (-3)*x[3] + (-3)*x[4] + (4)*x[5] + (-4)*x[6] + (5)*x[7];
s.t. c1: (x[2] >= 0) <==> (x[3] + x[2] <= 1);
s.t. c2: atleast 2 {i in 1..7} (x[i] >= 1);
s.t. c3: (if x[7] <= 1 then x[1] else x[6]) < 1;
s.t. c4: atmost 3 {i in 1..7} (x[i] >= 1);
