var x{1..5} integer >= -2 <= 0;

maximize o: -x[1] + 5*x[2] + 4*x[3] - 2*x[4] + 0*x[5];

s.t. c1: ((x[5] <= 0) ==> (x[5] > -2)) and ((x[5] > -2) ==> (x[3] >= 0));
s.t. c2: not ((x[2] != -2) and (x[4] >= 0));
s.t. c3: ((if x[1] >= -1 then 1 else 0) + (if x[2] >= -1 then 1 else 0)) >= 1;
s.t. c4: alldiff(x[2], x[3]);
