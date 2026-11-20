var x{1..5} integer >= 0 <= 3;

maximize o: 17*x[1] - 333*x[2] + 1000*x[3] - 1000*x[4] - 1000*x[5];

s.t. c1: (if x[2] != 3 then x[3] else x[3]) != 1;
s.t. c2: ((if x[4] <= 3 then 1 else 0) + (if x[1] >= 3 then 1 else 0)) >= 1;
s.t. c3: alldiff(x[4], x[5], x[2]);
s.t. c4: (if x[4] > 3 then x[4] else x[2]) < 3;
