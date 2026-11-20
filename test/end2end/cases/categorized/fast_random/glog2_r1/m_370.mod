var x{1..4} integer >= 0 <= 3;

maximize o: 2*x[1] + 2*x[2] - x[3] + x[4];

s.t. c1: alldiff(x[2], x[4]);
s.t. c2: (x[1] != 2) or (x[3] == 3) or (x[1] >= 0);
