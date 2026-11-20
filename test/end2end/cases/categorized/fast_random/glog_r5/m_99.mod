var x{1..6} integer >= 0 <= 3;
minimize o: (-3)*x[1] + (-3)*x[2] + (4)*x[3] + (-4)*x[4] + (4)*x[5] + (3)*x[6];
s.t. c1: not ((x[5] < 1) and (x[3] != 1));
s.t. c2: atmost 1 {i in 1..6} (x[i] >= 3);
s.t. c3: (if x[1] < 0 then x[1] else x[6]) < 2;
