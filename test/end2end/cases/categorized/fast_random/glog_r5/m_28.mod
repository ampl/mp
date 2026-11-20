var x{1..5} integer >= 0 <= 4;

maximize o: x[1] + 3*x[2] + x[3] + x[4] - 3*x[5];

s.t. c1: (if x[5] < 2 then x[2] else x[5]) > 2;
s.t. c2: (if x[2] > 4 then x[2] else x[5]) <= 2;
