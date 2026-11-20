var x{1..6} integer >= 0 <= 3;

minimize o: 0*x[1] + 0*x[2] + 0*x[3] + 0*x[4] + 0*x[5] + 0*x[6];

s.t. c1: ((if x[6] == 3 then 1 else 0) + (if x[4] >= 3 then 1 else 0)) >= 1;
s.t. c2: ((if x[3] > 3 then 1 else 0) + (if x[2] >= 3 then 1 else 0)) >= 1;
s.t. c3: max(x[6], x[3]) - min(x[2], x[6]) == 0;
s.t. c4: alldiff(x[3], x[6], x[5]);
