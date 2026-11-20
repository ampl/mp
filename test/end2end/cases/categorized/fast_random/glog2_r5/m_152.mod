var x{1..6} integer >= 0 <= 3;

maximize o: 2*x[1] + x[2] - 3*x[3] + 0*x[4] + 0*x[5] - x[6];

s.t. c1: atleast 3 {i in 1..6} (x[i] >= 0);
s.t. c2: not ((x[5] == 3) and (x[5] == 1));
s.t. c3: (x[1] <= 0) or (x[4] < 0) or (x[2] <= 1);
s.t. c4: alldiff(x[3], x[6]);
