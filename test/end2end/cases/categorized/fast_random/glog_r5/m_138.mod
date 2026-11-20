var x{1..5} integer >= -2 <= 1;
minimize o: (-5)*x[1] + (-3)*x[2] + (0)*x[3] + (4)*x[4] + (1)*x[5];
s.t. c1: (if x[3] <= 1 then x[3] else x[3]) > -1;
s.t. c2: atmost 2 {i in 1..5} (x[i] >= -1);
s.t. c3: (x[5] > -1) or (x[5] <= 1) or (x[1] != -2);
s.t. c4: (x[3] <= 0) or (x[4] != 1) or (x[5] > -2);
