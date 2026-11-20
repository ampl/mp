var x{1..3} integer >= 0 <= 4;

maximize o: -5*x[1] + 2*x[2] + 6*x[3];

s.t. c1: exists {i in 1..3} x[i] == 0;
s.t. c2: alldiff(x[2], x[1], x[3]);
s.t. c3: (x[3] > 0) and (x[1] > 3) ==> (x[3] >= 0);
