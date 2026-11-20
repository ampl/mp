var x{1..6} integer >= -2 <= 0;

minimize o: -5*x[1] + 3*x[2] - 4*x[3] - 5*x[4] + 6*x[5] + 0*x[6];

s.t. c1: numberof -3 in (x[1], x[2], x[3], x[4], x[5], x[6]) >= 3;
s.t. c2: alldiff(x[3], x[1], x[4]);
s.t. c3: (if x[3] != -3 then x[4] else x[1]) < 0;
s.t. c4: (x[2] > 0) ==> (x[1] >= -2) else (x[4] > 0);
