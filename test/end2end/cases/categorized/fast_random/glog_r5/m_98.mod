var x{1..4} integer >= 0 <= 3;

maximize o: 0*x[1] + 0*x[2] + 2*x[3] - 5*x[4];

s.t. c1: (x[2] == 1) ==> (x[3] <= 3);
s.t. c2: exists {i in 1..4} x[i] == 2;
s.t. c3: (x[2] <= 2) <==> (x[2] + x[4] >= 3);
s.t. c4: exists {i in 1..4} x[i] == 2;
