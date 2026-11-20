var x{1..6} integer >= -2 <= 1;

maximize o: 4*x[1] + x[2] - 5*x[3] + 0*x[4] - 4*x[5] + 2*x[6];

s.t. c1: not ((x[3] == 0) and (x[5] > -1));
s.t. c2: (if x[6] < -2 then x[1] else x[1]) == 0;
