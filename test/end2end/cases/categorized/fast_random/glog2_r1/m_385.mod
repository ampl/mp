var x{1..3} integer >= 0 <= 3;

maximize o: 4*x[1] + x[2] + 5*x[3];

s.t. c1: (if x[3] <= 2 then x[3] else x[2]) > 3;
s.t. c2: not ((x[1] >= 2) and (x[2] < 2));
s.t. c3: ((if x[2] != 0 then 1 else 0) + (if x[3] >= 0 then 1 else 0)) >= 1;
s.t. c4: alldiff(x[2], x[3], x[1]);
s.t. c5: abs(x[2] - x[2]) > 0;
s.t. c6: (x[3] < 3) <==> (x[2] + x[3] >= 0);
