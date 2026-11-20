var x{1..4} integer >= 0 <= 2;

minimize o: 3000*x[1] + 17*x[2] + 34*x[3] - 666*x[4];

s.t. c1: exists {i in 1..4} x[i] == 1;
