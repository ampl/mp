var x{1..5} integer >= -2 <= 0;
minimize o: (0)*x[1] + (0)*x[2] + (0)*x[3] + (0)*x[4] + (0)*x[5];
s.t. c1: (if x[4] != -1 then x[5] else x[3]) >= -2;
s.t. c2: (if x[5] > -2 then x[3] else x[2]) < 0;
