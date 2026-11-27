var x{1..6} integer >= 0 <= 4;

minimize o: 51*x[1] - 3000*x[2] - 666*x[3] + 34*x[4] + 51*x[5] - 999*x[6];

s.t. c1: exists {i in 1..6} x[i] == 2;
