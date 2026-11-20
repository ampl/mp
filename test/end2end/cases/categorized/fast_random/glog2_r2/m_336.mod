var x{1..6} integer >= -2 <= 1;
minimize o: (0)*x[1] + (0)*x[2] + (0)*x[3] + (0)*x[4] + (0)*x[5] + (0)*x[6];
s.t. c1: (if x[1] >= -2 then x[4] else x[1]) > -2;
