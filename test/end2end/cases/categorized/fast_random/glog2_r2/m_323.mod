var x{1..5} integer >= 0 <= 4;

maximize o: 0*x[1] + 0*x[2] - 2*x[3] + 4*x[4] - 2*x[5];

s.t. c1: (x[1] <= 4) <==> (x[1] + x[1] != 1);
s.t. c2: not ((x[2] == 3) and (x[4] < 0));
s.t. c3: not ((x[3] >= 3) and (x[1] == 1));
s.t. c4: exists {i in 1..5} x[i] == 3;
s.t. c5: x[4] == 4;
