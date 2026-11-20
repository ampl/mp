var x{1..5} integer >= -2 <= 2;

maximize o: 0*x[1] + 4*x[2] - 5*x[3] + 0*x[4] + 5*x[5];

s.t. c1: (x[4] <= -1) <==> (x[5] + x[3] < 0);
s.t. c2: max(x[4], x[2]) - min(x[4], x[1]) == -1;
s.t. c3: exists {i in 1..5} x[i] == -1;
