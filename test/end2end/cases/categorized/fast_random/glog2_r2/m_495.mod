var x{1..4} integer >= -2 <= 2;

maximize o: -5*x[1] - 2*x[2] - 4*x[3] + 4*x[4];

s.t. c1: exists {i in 1..4} x[i] == -1;
s.t. c2: atleast 2 {i in 1..4} (x[i] >= 0);
s.t. c3: alldiff(x[1], x[2], x[3]);
s.t. c4: alldiff(x[3], x[4]);
s.t. c5: (x[2] < -1) <==> (x[4] + x[4] >= -1);
