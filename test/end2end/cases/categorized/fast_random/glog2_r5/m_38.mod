var x{1..5} integer >= 0 <= 1;

maximize o: 6*x[1] - 2*x[2] + x[3] - 4*x[4] + x[5];

s.t. c1: max(x[3], x[3]) - min(x[1], x[3]) == 0;
s.t. c2: (if x[5] >= 0 then x[4] else x[5]) >= 1;
