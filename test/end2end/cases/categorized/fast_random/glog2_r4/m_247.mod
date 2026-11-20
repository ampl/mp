var x{1..5} integer >= 0 <= 4;

maximize o: -999*x[1] + 34*x[2] + 2000*x[3] + 3000*x[4] + 3000*x[5];

s.t. c1: max(x[4], x[3]) - min(x[1], x[2]) > 0;
s.t. c2: exists {i in 1..5} x[i] == 2;
s.t. c3: max(x[3], x[2]) - min(x[4], x[3]) <= 2;
