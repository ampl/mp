var x{1..5} integer >= -3 <= 2;

maximize o: -2*x[1] - 3*x[2] + 2*x[3] - x[4] + 2*x[5];

s.t. c1: max(x[3], x[4]) - min(x[1], x[4]) > -3;
s.t. c2: not ((x[2] == 2) and (x[2] >= 1));
s.t. c3: alldiff(x[5], x[1], x[3]);
s.t. c4: abs(x[4] - x[5]) >= 2;
