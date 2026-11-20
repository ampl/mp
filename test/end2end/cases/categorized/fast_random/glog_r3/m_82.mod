var x{1..5} integer >= 0 <= 4;

minimize o: 3*x[1] - 4*x[2] + 0*x[3] + 3*x[4] - x[5];

s.t. c1: exists {i in 1..5} x[i] == 2;
s.t. c2: ((x[5] < 4) ==> (x[4] <= 1)) and ((x[4] <= 1) ==> (x[4] >= 4));
