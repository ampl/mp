var x{1..5} integer >= 0 <= 4;

minimize o: 34*x[1] + 34*x[2] - 1000*x[3] + 34*x[4] + 17*x[5];

s.t. c1: max(x[4], x[2]) - min(x[4], x[3]) > 2;
s.t. c2: exists {i in 1..5} x[i] == 4;
