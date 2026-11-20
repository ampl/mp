var x{1..5} integer >= 0 <= 3;

minimize o: 34*x[1] + 1000*x[2] - 2000*x[3] - 999*x[4] + 1000*x[5];

s.t. c1: (x[3] < 3) or (x[5] > 2) or (x[1] != 1);
s.t. c2: max(x[2], x[3]) - min(x[4], x[3]) > 3;
s.t. c3: alldiff(x[3], x[4], x[2]);
