var x{1..6} integer >= -2 <= 1;

maximize o: -2*x[1] + 4*x[2] + 3*x[3] - 5*x[4] - 5*x[5] + 0*x[6];

s.t. c1: ((if x[3] <= -2 then 1 else 0) + (if x[3] >= -2 then 1 else 0)) >= 1;
s.t. c2: (x[2] < 1) or (x[6] > 1) or (x[4] >= -1);
s.t. c3: alldiff(x[4], x[3], x[5], x[6]);
s.t. c4: ((if x[3] <= -2 then 1 else 0) + (if x[3] >= -2 then 1 else 0)) >= 1;
s.t. c5: alldiff(x[4], x[3], x[5], x[6]);
