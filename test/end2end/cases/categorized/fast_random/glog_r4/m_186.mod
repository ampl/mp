var x{1..4} integer >= -2 <= 1;

maximize o: x[1] + x[2] - x[3] - 3*x[4];

s.t. c1: (x[2] <= -2) <==> (x[3] + x[2] <= -1);
s.t. c2: exists {i in 1..4} x[i] == -1;
s.t. c3: exists {i in 1..4} x[i] == 0;
s.t. c4: max(x[3], x[3]) - min(x[4], x[2]) != 0;
