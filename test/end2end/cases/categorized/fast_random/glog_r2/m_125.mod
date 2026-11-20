var x{1..5} integer >= 0 <= 3;

minimize o: -x[1] - 4*x[2] + 4*x[3] - 2*x[4] - 2*x[5];

s.t. c1: (x[1] <= 0) <==> (x[2] + x[3] >= 2);
s.t. c2: (x[1] > 1) or (x[3] != 1) or (x[3] >= 3);
s.t. c3: exists {i in 1..5} x[i] == 0;
