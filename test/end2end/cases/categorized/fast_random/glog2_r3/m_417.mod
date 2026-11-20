var x{1..6} integer >= 0 <= 3;

maximize o: 34*x[1] + 17*x[2] + 1000*x[3] - 3000*x[4] - 2000*x[5] - 3000*x[6];

s.t. c1: numberof 2 in (x[1], x[2], x[3], x[4], x[5], x[6]) < 2;
s.t. c2: alldiff(x[3], x[1], x[4], x[6]);
s.t. c3: max(x[2], x[1]) - min(x[2], x[1]) > 3;
s.t. c4: exists {i in 1..6} x[i] == 2;
