var x{1..3} integer >= 0 <= 4;

maximize o: -4*x[1] + 6*x[2] + 6*x[3];

s.t. c1: (x[2] != 0) ==> (x[2] == 3) else (x[3] == 0);
s.t. c2: (x[1] <= 1) <==> (x[2] + x[2] <= 1);
s.t. c3: alldiff(x[2], x[1]);
s.t. c4: (x[2] != 0) ==> (x[2] == 3) else (x[3] == 0);
s.t. c5: alldiff(x[2], x[1]);
