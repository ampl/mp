var x{1..3} integer >= 0 <= 3;

maximize o: 6*x[1] + 3*x[2] - 4*x[3];

s.t. c1: ((if x[1] >= 3 then 1 else 0) + (if x[3] >= 3 then 1 else 0)) >= 1;
s.t. c2: alldiff(x[3], x[1], x[2]);
s.t. c3: not ((x[1] == 3) and (x[1] >= 0));
s.t. c4: max(x[2], x[1]) - min(x[2], x[3]) != 0;
s.t. c5: (x[2] <= 3) or (x[2] < 1) or (x[3] <= 0);
