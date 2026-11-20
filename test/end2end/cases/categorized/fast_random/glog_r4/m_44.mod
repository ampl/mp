var x{1..6} integer >= 0 <= 2;
minimize o: (1)*x[1] + (4)*x[2] + (5)*x[3] + (4)*x[4] + (1)*x[5] + (-2)*x[6];
s.t. c1: (if x[5] <= 2 then x[6] else x[6]) < 1;
