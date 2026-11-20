var x{1..6} integer >= 0 <= 3;

maximize o: 4*x[1] - 4*x[2] - 3*x[3] - 4*x[4] - 3*x[5] + 5*x[6];

s.t. c1: max(x[6], x[4]) - min(x[3], x[6]) >= 3;
s.t. c2: (x[4] < 2) or (x[1] < 3) or (x[3] != 2);
s.t. c3: alldiff(x[1], x[2], x[6]);
s.t. c4: max(x[6], x[4]) - min(x[3], x[6]) >= 3;
s.t. c5: alldiff(x[1], x[2], x[6]);
