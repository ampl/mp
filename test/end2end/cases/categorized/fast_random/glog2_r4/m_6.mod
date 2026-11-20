var x{1..3} integer >= 0 <= 4;

maximize o: 0*x[1] + 0*x[2] + 0*x[3];

s.t. c1: not ((x[2] != 2) and (x[2] != 1));
s.t. c2: alldiff(x[2], x[3], x[1]);
s.t. c3: ((if x[2] < 4 then 1 else 0) + (if x[2] >= 4 then 1 else 0)) >= 1;
