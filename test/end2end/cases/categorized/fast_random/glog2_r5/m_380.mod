var x{1..5} integer >= -2 <= 1;

maximize o: -3*x[1] + x[2] + 3*x[3] + 0*x[4] + 2*x[5];

s.t. c1: alldiff(x[1], x[4], x[3]);
s.t. c2: (if x[3] < -2 then x[1] else x[1]) != -1;
s.t. c3: (x[1] <= 0) or (x[1] >= 0) or (x[2] <= 1);
s.t. c4: ((if x[4] <= 1 then 1 else 0) + (if x[3] >= 1 then 1 else 0)) >= 1;
