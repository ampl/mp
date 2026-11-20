var x{1..5} integer >= 0 <= 3;

maximize o: -4*x[1] + 0*x[2] - 4*x[3] - 2*x[4] + 3*x[5];

s.t. c1: ((if x[1] != 2 then 1 else 0) + (if x[3] >= 2 then 1 else 0)) >= 1;
s.t. c2: (if x[4] == 1 then x[5] else x[4]) >= 3;
s.t. c3: (if x[3] == 2 then x[1] else x[5]) < 3;
