var x{1..4} integer >= -2 <= 1;

maximize o: -5*x[1] + 0*x[2] - 3*x[3] + 3*x[4];

s.t. c1: alldiff(x[1], x[2], x[4], x[3]);
s.t. c2: atleast 1 {i in 1..4} (x[i] >= -1);
s.t. c3: (x[3] != 1) ==> (x[3] < -1) else (x[3] > 1);
