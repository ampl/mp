var x{1..5} integer >= 0 <= 4;

maximize o: -2*x[1] - 3*x[2] + 0*x[3] + 3*x[4] + 0*x[5];

s.t. c1: (if x[2] >= 4 then x[5] else x[5]) > -1;
s.t. c2: not ((x[4] <= 5) and (x[5] > 5));
s.t. c3: not ((x[4] > 0) and (x[1] >= -1));
