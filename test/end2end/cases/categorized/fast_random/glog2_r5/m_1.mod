var x{1..4} integer >= 0 <= 4;

maximize o: x[1] - 4*x[2] - x[3] + 4*x[4];

s.t. c1: exists {i in 1..4} x[i] == -1;
s.t. c2: max(x[2], x[1]) - min(x[2], x[1]) < 5;
s.t. c3: exists {i in 1..4} x[i] == 0;
s.t. c4: alldiff(x[2], x[1]);
