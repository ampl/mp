var x{1..6} integer >= 0 <= 4;

maximize o: -4*x[1] + 4*x[2] + 0*x[3] + 0*x[4] - 5*x[5] + 6*x[6];

s.t. c1: max(x[3], x[5]) - min(x[4], x[4]) != 1;
s.t. c2: (if x[5] != 2 then x[5] else x[4]) == 3;
s.t. c3: x[1] == 4;
