var x{1..6} integer >= 0 <= 2;

maximize o: -4*x[1] - 3*x[2] - 4*x[3] - 4*x[4] + 3*x[5] + 2*x[6];

s.t. c1: (if x[6] != 1 then x[2] else x[1]) != 1;
s.t. c2: x[6] == 2;
