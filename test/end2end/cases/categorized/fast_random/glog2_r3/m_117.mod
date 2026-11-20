var x{1..6} integer >= 0 <= 2;

maximize o:
   -999*x[1] + 1000*x[2] - 2000*x[3] - 2000*x[4] + 3000*x[5] - 2000*x[6];

s.t. c1: (if x[1] <= 1 then x[5] else x[5]) <= 2;
s.t. c2: (x[6] <= 2) and (x[4] < 0) ==> (x[3] >= 2);
s.t. c3: ((if x[4] > 0 then 1 else 0) + (if x[1] >= 0 then 1 else 0)) >= 1;
s.t. c4: alldiff(x[1], x[5], x[6]);
