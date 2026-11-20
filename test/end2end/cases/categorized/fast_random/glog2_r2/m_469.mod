var x{1..3} integer >= 0 <= 4;

minimize o: 6*x[1] + 5*x[2] + x[3];

s.t. c1: alldiff(x[1], x[2]);
s.t. c2: ((if x[1] >= 2 then 1 else 0) + (if x[2] >= 2 then 1 else 0)) >= 1;
s.t. c3: alldiff(x[1], x[2]);
s.t. c4: ((if x[1] >= 2 then 1 else 0) + (if x[2] >= 2 then 1 else 0)) >= 1;
