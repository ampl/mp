var x{1..4} integer >= -2 <= 2;

maximize o: -5*x[1] + 3*x[2] + 2*x[3] - 3*x[4];

s.t. c1: exists {i in 1..4} x[i] == 1;
s.t. c2: max(x[3], x[3]) - min(x[2], x[2]) >= 1;
s.t. c3: alldiff(x[3], x[1]);
s.t. c4: (x[3] >= 0) ==> (x[3] != 0);
