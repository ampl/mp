var x{1..5} integer >= -2 <= 1;

minimize o: -2000*x[1] + 17*x[2] - 666*x[3] + 3000*x[4] - 1000*x[5];

s.t. c1: exists {i in 1..5} x[i] == -2;
s.t. c2: max(x[4], x[4]) - min(x[3], x[1]) > -2;
s.t. c3: exists {i in 1..5} x[i] == -2;
