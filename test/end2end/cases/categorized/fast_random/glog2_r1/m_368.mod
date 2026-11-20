var x{1..6} integer >= 0 <= 1;
minimize o: (1)*x[1] + (-1)*x[2] + (-5)*x[3] + (4)*x[4] + (-5)*x[5] + (-2)*x[6];
s.t. c1: (if x[2] > 1 then x[5] else x[4]) >= 1;
s.t. c2: (x[2] == 0) or (x[6] != 1) or (x[5] == 0);
