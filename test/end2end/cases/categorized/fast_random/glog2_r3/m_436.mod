var x{1..5} integer >= -2 <= 1;

maximize o: 0*x[1] + 0*x[2] + 0*x[3] + 0*x[4] + 0*x[5];

s.t. c1: abs(x[1] - x[4]) >= 0;
s.t. c2: ((if x[2] > -1 then 1 else 0) + (if x[5] >= -1 then 1 else 0)) >= 1;
s.t. c3: not ((x[2] < 1) and (x[1] == 1));
