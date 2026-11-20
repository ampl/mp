var x{1..5} integer >= 0 <= 3;

maximize o: -2000*x[1] + 2000*x[2] + 34*x[3] - 1000*x[4] - 666*x[5];

s.t. c1: max(x[1], x[4]) - min(x[1], x[5]) == 0;
s.t. c2: max(x[1], x[2]) - min(x[4], x[2]) <= 3;
