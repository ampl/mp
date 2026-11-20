var x{1..3} integer >= 0 <= 3;

maximize o: -x[1] - 5*x[2] + 6*x[3];

s.t. c1: (if x[3] > 4 then x[3] else x[3]) <= 0;
s.t. c2: (x[1] != -1) or (x[1] >= -1) or (x[3] > 0);
