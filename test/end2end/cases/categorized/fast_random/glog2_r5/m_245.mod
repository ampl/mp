var x{1..5} integer >= 0 <= 4;

minimize o: 0*x[1] - 3*x[2] + 3*x[3] + 0*x[4] + 3*x[5];

s.t. c1: exists {i in 1..5} x[i] == 1;
s.t. c2: (x[5] <= 2) and (x[3] <= 4) ==> (x[5] >= 2);
s.t. c3: not ((x[1] >= 2) and (x[2] >= 4));
s.t. c4: abs(x[1] - x[5]) < 1;
s.t. c5: (x[4] < 1) ==> (x[2] > 1);
