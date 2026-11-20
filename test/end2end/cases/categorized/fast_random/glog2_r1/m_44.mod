var x{1..3} integer >= -3 <= 2;

minimize o: -5*x[1] + x[2] + 3*x[3];

s.t. c1: not ((x[3] >= 0) and (x[3] >= -1));
s.t. c2: abs(x[3] - x[3]) >= 2;
s.t. c3: alldiff(x[3], x[1], x[2]);
s.t. c4: ((if x[3] < -2 then 1 else 0) + (if x[3] >= -2 then 1 else 0)) >= 1;
