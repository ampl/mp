var x{1..5} integer >= -2 <= 0;

maximize o: 17*x[1] - 3000*x[2] + 34*x[3] + 2000*x[4] - 1000*x[5];

s.t. c1: max(x[2], x[4]) - min(x[1], x[5]) == 0;
s.t. c2: not ((x[3] <= -1) and (x[1] <= -2));
s.t. c3: not ((x[5] >= -2) and (x[2] > -2));
s.t. c4: abs(x[4] - x[3]) > 1;
