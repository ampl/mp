var x{1..5} integer >= 0 <= 2;

maximize o: -x[1] - 4*x[2] + 6*x[3] - 5*x[4] - 5*x[5];

s.t. c1: alldiff(x[1], x[3]);
s.t. c2: (x[4] > 1) <==> (x[2] + x[3] != 2);
s.t. c3: (x[5] > 0) ==> (x[5] < 2) else (x[5] == 1);
s.t. c4: not ((x[4] < 2) and (x[4] == 2));
