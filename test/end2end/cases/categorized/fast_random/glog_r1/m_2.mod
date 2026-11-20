var x{1..4} integer >= 0 <= 4;

maximize o: -2*x[1] + x[2] - 2*x[3] - x[4];

s.t. c1: (x[3] <= 4) ==> (x[2] < 1) else (x[3] > 2);
s.t. c2: alldiff(x[1], x[4], x[2], x[3]);
s.t. c3: alldiff(x[3], x[2]);
