var x{1..5} integer >= 0 <= 4;

minimize o: -5*x[1] + 4*x[2] - x[3] - 5*x[4] + 2*x[5];

s.t. c1: abs(x[4] - x[1]) > 2;
s.t. c2: exists {i in 1..5} x[i] == 1;
s.t. c3: (x[4] >= 4) <==> (x[2] + x[5] <= 3);
