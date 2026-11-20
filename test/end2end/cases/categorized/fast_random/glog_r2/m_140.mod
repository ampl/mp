var x{1..3} integer >= 0 <= 3;

maximize o: -2*x[1] + 3*x[2] - x[3];

s.t. c1: max(x[1], x[3]) - min(x[2], x[3]) < 3;
s.t. c2: (x[1] < 3) or (x[3] < 0) or (x[3] != 0);
s.t. c3: alldiff(x[3], x[2]);
s.t. c4: exists {i in 1..3} x[i] == 2;
