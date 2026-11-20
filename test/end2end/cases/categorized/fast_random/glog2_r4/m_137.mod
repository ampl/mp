var x{1..6} integer >= 0 <= 3;

minimize o: -333*x[1] + 3000*x[2] - 666*x[3] - 1000*x[4] - 333*x[5] + 51*x[6];

s.t. c1: exists {i in 1..6} x[i] == 1;
s.t. c2: (x[2] > 0) ==> (x[2] < 1) else (x[5] != 0);
