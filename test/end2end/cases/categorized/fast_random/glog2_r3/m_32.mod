var x{1..6} integer >= -2 <= 2;

maximize o: x[1] - 3*x[2] - 2*x[3] + 6*x[4] + 4*x[5] - 5*x[6];

s.t. c1: not ((x[6] >= 0) and (x[1] < 2));
s.t. c2: (if x[1] <= 1 then x[1] else x[1]) > 2;
