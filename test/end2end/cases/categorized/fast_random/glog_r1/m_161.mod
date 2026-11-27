var x{1..3} integer >= 0 <= 3;

maximize o: 0*x[1] + 0*x[2] - 3*x[3];

s.t. c1: alldiff(x[3], x[1]);
s.t. c2: exists {i in 1..3} x[i] == 1;
s.t. c3: (x[1] == 2) ==> (x[3] <= 2) else (x[3] <= 0);
