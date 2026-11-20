var x{1..6} integer >= -2 <= 1;

maximize o: 0*x[1] + 6*x[2] + 0*x[3] - 4*x[4] - x[5] + x[6];

s.t. c1: (if x[3] != -2 then x[3] else x[2]) >= -1;
