var x{1..6} integer >= 0 <= 2;

minimize o: -5*x[1] - 4*x[2] - 5*x[3] + x[4] + 4*x[5] + 5*x[6];

s.t. c1: ((x[4] >= 1) ==> (x[1] < 0)) and ((x[1] < 0) ==> (x[6] >= 1));
s.t. c2: exists {i in 1..6} x[i] == 1;
s.t. c3: abs(x[2] - x[6]) > 2;
