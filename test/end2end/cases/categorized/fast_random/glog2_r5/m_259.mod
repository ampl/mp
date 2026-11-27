var x{1..4} integer >= 0 <= 2;

maximize o: x[1] + 6*x[2] + 4*x[3] + 5*x[4];

s.t. c1: not ((x[3] < 2) and (x[1] < 0));
s.t. c2: (x[1] <= 2) ==> (x[2] <= 0);
s.t. c3: alldiff(x[4], x[1], x[3]);
s.t. c4: alldiff(x[4], x[2], x[3]);
s.t. c5: not ((x[3] < 2) and (x[1] < 0));
s.t. c6: alldiff(x[4], x[2], x[3]);
