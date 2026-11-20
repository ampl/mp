var x{1..6} integer >= -2 <= 2;

maximize o: x[1] - 4*x[2] - 2*x[3] + 2*x[4] + 4*x[5] + x[6];

s.t. c1: exists {i in 1..6} x[i] == -2;
s.t. c2: abs(x[2] - x[1]) == 4;
s.t. c3: (x[6] != 0) ==> (x[2] > 2);
s.t. c4: (x[1] >= -1) ==> (x[1] < -2);
