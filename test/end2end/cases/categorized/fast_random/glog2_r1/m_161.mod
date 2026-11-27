var x{1..6} integer >= 0 <= 3;

maximize o: -x[1] - 5*x[2] + x[3] - 5*x[4] - 3*x[5] + 3*x[6];

s.t. c1: alldiff(x[5], x[1], x[4]);
s.t. c2: (x[3] >= 3) or (x[1] > 4) or (x[5] <= 3);
s.t. c3: not ((x[5] > 4) and (x[5] == -1));
s.t. c4: max(x[4], x[5]) - min(x[3], x[2]) > -1;
