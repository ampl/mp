var x{1..5} integer >= 0 <= 2;

maximize o: -2*x[1] - 2*x[2] - 4*x[3] + x[4] + 2*x[5];

s.t. c1: (if x[4] == 2 then x[1] else x[5]) < 2;
s.t. c2: (if x[1] != 0 then x[5] else x[5]) <= 2;
