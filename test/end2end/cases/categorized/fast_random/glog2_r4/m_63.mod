var x{1..6} integer >= 0 <= 2;

maximize o: -3*x[1] + 3*x[2] - 2*x[3] - 5*x[4] + 0*x[5] + x[6];

s.t. c1: (x[5] < 0) <==> (x[2] + x[3] == 0);
s.t. c2: alldiff(x[2], x[6], x[4]);
s.t. c3: ((if x[6] == 1 then 1 else 0) + (if x[4] >= 1 then 1 else 0)) >= 1;
s.t. c4: (x[2] != 0) <==> (x[1] + x[6] <= 2);
s.t. c5: x[6] == 0;
