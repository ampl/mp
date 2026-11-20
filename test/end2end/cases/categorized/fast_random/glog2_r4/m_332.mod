var x{1..3} integer >= 0 <= 2;

maximize o: 0*x[1] - x[2] + 4*x[3];

s.t. c1: alldiff(x[1], x[2], x[3]);
s.t. c2: exists {i in 1..3} x[i] == 2;
s.t. c3: (x[2] > 2) or (x[3] >= 2) or (x[2] > 1);
s.t. c4: (x[3] < 2) ==> (x[3] < 0) else (x[3] == 1);
