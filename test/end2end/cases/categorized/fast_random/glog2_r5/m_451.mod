var x{1..6} integer >= 0 <= 4;

maximize o: -5*x[1] + x[2] + 0*x[3] + x[4] - x[5] + 0*x[6];

s.t. c1: (x[3] != 5) or (x[4] >= -1) or (x[2] >= -1);
s.t. c2: exists {i in 1..6} x[i] == -1;
s.t. c3: (x[5] == 5) <==> (x[5] + x[4] >= 0);
