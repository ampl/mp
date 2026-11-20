var x{1..5} integer >= -2 <= 2;

maximize o: 2*x[1] + 0*x[2] - 2*x[3] - 2*x[4] - 4*x[5];

s.t. c1: numberof 2 in (x[1], x[2], x[3], x[4], x[5]) > 0;
s.t. c2: alldiff(x[5], x[2], x[1], x[4]);
s.t. c3: ((if x[5] < -1 then 1 else 0) + (if x[5] >= -1 then 1 else 0)) >= 1;
s.t. c4: (x[4] <= -2) ==> (x[3] == 1) else (x[1] < -2);
