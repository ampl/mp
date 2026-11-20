var x{1..6} integer >= 0 <= 3;

maximize o: -2*x[1] + 0*x[2] + 6*x[3] + 6*x[4] - 4*x[5] + 6*x[6];

s.t. c1: (if x[3] < 1 then x[6] else x[3]) > 0;
