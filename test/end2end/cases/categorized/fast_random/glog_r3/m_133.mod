var x{1..4} integer >= -2 <= 0;

minimize o: -5*x[1] + 6*x[2] - 4*x[3] - 3*x[4];

s.t. c1: alldiff(x[1], x[3]);
s.t. c2: ((if x[2] < -1 then 1 else 0) + (if x[1] >= -1 then 1 else 0)) >= 1;
s.t. c3: (if x[2] < -2 then x[2] else x[4]) > -2;
s.t. c4: alldiff(x[4], x[1], x[2], x[3]);
