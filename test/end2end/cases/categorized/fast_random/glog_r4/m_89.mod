var x{1..4} integer >= -2 <= 0;

maximize o: -3*x[1] + 5*x[2] - x[3] + 4*x[4];

s.t. c1: (x[2] <= -1) <==> (x[4] + x[1] <= 0);
s.t. c2: exists {i in 1..4} x[i] == 0;
s.t. c3: (x[4] <= -2) ==> (x[4] >= 0) else (x[2] <= 0);
s.t. c4: exists {i in 1..4} x[i] == 0;
