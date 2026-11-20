var x{1..4} integer >= 0 <= 4;
minimize o: (-4)*x[1] + (-5)*x[2] + (4)*x[3] + (4)*x[4];
s.t. c1: exists {i in 1..4} x[i] == 0;
s.t. c2: (if x[2] < 0 then x[4] else x[2]) <= 1;
