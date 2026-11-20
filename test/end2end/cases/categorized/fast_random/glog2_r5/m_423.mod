var x{1..4} integer >= -2 <= 0;

maximize o: 3*x[1] + 4*x[2] + 4*x[3] - 5*x[4];

s.t. c1: ((if x[1] == 0 then 1 else 0) + (if x[4] >= 0 then 1 else 0)) >= 1;
s.t. c2: x[3] == -2;
