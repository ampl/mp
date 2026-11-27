var x{1..5} integer >= 0 <= 3;

minimize o: 6*x[1] + 3*x[2] - 5*x[3] - 4*x[4] + 6*x[5];

s.t. c1: ((if x[5] > 3 then 1 else 0) + (if x[2] >= 3 then 1 else 0)) >= 1;
s.t. c2: numberof 2 in (x[1], x[2], x[3], x[4], x[5]) < 1;
s.t. c3: alldiff(x[2], x[3], x[4], x[1]);
s.t. c4: ((if x[5] > 3 then 1 else 0) + (if x[2] >= 3 then 1 else 0)) >= 1;
s.t. c5: alldiff(x[2], x[3], x[4], x[1]);
