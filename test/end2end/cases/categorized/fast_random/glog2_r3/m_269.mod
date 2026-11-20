var x{1..4} integer >= 0 <= 4;

maximize o: 4*x[1] - 3*x[2] - 4*x[3] + x[4];

s.t. c1: exists {i in 1..4} x[i] == 1;
s.t. c2: not ((x[1] < 0) and (x[3] >= 1));
s.t. c3: (x[4] <= 4) <==> (x[1] + x[2] >= 0);
s.t. c4: exists {i in 1..4} x[i] == 2;
s.t. c5: exists {i in 1..4} x[i] == 1;
s.t. c6: exists {i in 1..4} x[i] == 2;
