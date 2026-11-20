var x{1..4} integer >= 0 <= 3;

minimize o: -3*x[1] + 3*x[2] + 2*x[3] + 4*x[4];

s.t. c1: alldiff(x[4], x[1], x[3], x[2]);
s.t. c2: max(x[3], x[2]) - min(x[2], x[4]) > 0;
s.t. c3: alldiff(x[4], x[1]);
s.t. c4: ((if x[2] != 3 then 1 else 0) + (if x[3] >= 3 then 1 else 0)) >= 1;
s.t. c5: alldiff(x[4], x[1], x[3], x[2]);
s.t. c6: ((if x[2] != 3 then 1 else 0) + (if x[3] >= 3 then 1 else 0)) >= 1;
