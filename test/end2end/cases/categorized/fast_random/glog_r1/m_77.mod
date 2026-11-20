var x{1..4} integer >= 0 <= 2;

maximize o: -2*x[1] + 5*x[2] - 2*x[3] + 2*x[4];

s.t. c1: ((if x[1] != 0 then 1 else 0) + (if x[3] >= 0 then 1 else 0)) >= 1;
s.t. c2: (x[3] == 0) or (x[4] != 2) or (x[2] != 0);
