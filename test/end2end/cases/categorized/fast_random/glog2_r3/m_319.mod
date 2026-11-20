var x{1..4} integer >= 0 <= 4;

maximize o: 6*x[1] - 2*x[2] + 2*x[3] - 4*x[4];

s.t. c1: alldiff(x[1], x[3]);
s.t. c2: atleast 2 {i in 1..4} (x[i] >= 2);
s.t. c3: alldiff(x[1], x[3]);
s.t. c4: atleast 2 {i in 1..4} (x[i] >= 2);
