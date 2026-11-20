var x{1..3} integer >= 0 <= 4;

minimize o: 3*x[1] + 4*x[2] - 2*x[3];

s.t. c1: not ((x[2] <= 0) and (x[1] >= -1));
s.t. c2: alldiff(x[1], x[3]);
s.t. c3: (x[3] != 5) <==> (x[1] + x[3] <= 5);
s.t. c4: (x[1] <= 5) ==> (x[2] > 0) else (x[1] < 4);
