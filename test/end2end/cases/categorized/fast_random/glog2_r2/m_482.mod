var x{1..3} integer >= 0 <= 3;

maximize o: 2*x[1] + 6*x[2] + 6*x[3];

s.t. c1: (x[1] >= 2) or (x[3] > 1) or (x[3] <= 0);
s.t. c2: alldiff(x[2], x[3], x[1]);
s.t. c3: (x[3] > 2) or (x[2] < 2) or (x[1] < 2);
